ModelTree_class_name <- "ModelTree"

generate.id <- function(existing.ids, prefix=NULL) {
  id <- NA
  counter <- 0
  while(is.na(id) | !is.na(match(id, existing.ids))) {
    counter <- counter + 1
    id <- paste(prefix, counter, sep=".")
  }
  return(id)
}

kern.sum <- function(a, b) {
  return(a + b)
}

kern.prod <- function(a, b) {
  return(a * b)
}

sum.grad <- function(a.grad, b.grad, ...) {
  return(a.grad + b.grad)
}

prod.grad <- function(a, b, a.grad, b.grad) {
  a1 <- b1 <- b.grad
  for (i in 1:dim(b.grad)[3]) {
    a1[, , i] <- a
    b1[, , i] <- b
  }
  return(a1 * b.grad + b1 * a.grad)
}

sum.hess <- function(a.hess, b.hess, ...) {
  return(a.hess + b.hess)
}

prod.hess <- function(a, b, a.grad, b.grad, a.hess, b.hess) {
  a1 <- b1 <- b.hess
  for (i in 1:dim(b.grad)[3]) {
    for (j in 1:dim(b.grad)[3]) {
        a1[, , i, j] <- a
        b1[, , i, j] <- b
    }
  }

  a.grad1 <- a.grad2 <- b.grad1 <- b.grad2 <- b.hess
  for (i in 1:dim(b.grad)[3]) {
    a.grad1[, , i, ] <- a.grad2[, , , i] <- a.grad
    b.grad1[, , i, ] <- b.grad2[, , , i] <- b.grad
  }

  return(a.hess * b1 + a.grad1 * b.grad2 + a.grad2 * b.grad1 + a1 * b.hess)
}

#' Create an Empty Model Tree
#'
#' Creates an instance of a modelTree object with no associated kernels.
#'
#' @return a modelTree object
#' @export
#'
#' @examples
#' mt <- create.model.tree()
create.model.tree <- function() {
  model.tree <- list()
  class(model.tree) <- ModelTree_class_name

  model.tree$operation.functions = list("+"=kern.sum,
                                        "*"=kern.prod
  )

  model.tree$operation.grad.functions = list("+"=sum.grad,
                                             "*"=prod.grad
  )

  model.tree$operation.hess.functions = list("+"=sum.hess,
                                             "*"=prod.hess
  )

  model.tree$tree <- data.frame(leftDaughter=numeric(0),
                                rightDaughter=numeric(0),
                                kernelInstanceID=character(0),
                                operationID=character(0),
                                stringsAsFactors=FALSE
  )

  model.tree$kernel.instances <- data.frame(kernelInstanceID=character(0),
                                            kernelClassName=character(0),
                                            stringsAsFactors=FALSE
  )

  model.tree$all.hyper.params <- c()

  # These are accessed by string IDs to make sure there's no issue
  # with accidentally accessing the nth element when we want the element
  # whose ID is "n".

  # This gives the indices for the hyperparameters in the hyperparameters
  # vector for a specific instance of a kernel. Indexed by kernelInstanceID.
  model.tree$kernel.inst.hyper.param.indices <- list()

  # This holds the base kernels (of class Kernel) available to the tree.
  # Indexed by kernelClassName.
  model.tree$kernel.objects <- list()

  return(model.tree)
}

find.root.node <- function(model.tree) {
  num.nodes <- nrow(model.tree$tree)
  if (num.nodes == 0) {
    root.node <- 1
  } else {
    child.nodes <- na.omit(c(model.tree$tree$leftDaughter, model.tree$tree$rightDaughter))
    root.node <- which(!1:num.nodes %in% child.nodes)
  }
  if (length(root.node) == 0){
    stop("No root node found.")
  }
  if (length(root.node) > 1){
    stop("Multiple root nodes found.")
  }
  return(root.node)
}

find.parent.node <- function(model.tree, node) {
  parent.node <- c(which(node == model.tree$tree$leftDaughter),
                   which(node == model.tree$tree$rightDaughter)
  )
  if (length(parent.node) > 1) {
    stop("Multiple parents found.")
  }
  if (length(parent.node) == 0) {
    # Presumably this is the root node. Return NULL.
    parent.node <- NULL
  }
  return(parent.node)
}

node.to.string <- function(model.tree, node, show.node.ids=FALSE) {
  if (nrow(model.tree$tree) == 0) {
    ret <- NULL
  } else {
    op <- model.tree$tree[node, "operationID"]
    ld <- model.tree$tree[node, "leftDaughter"]
    rd <- model.tree$tree[node, "rightDaughter"]
    if (is.na(op)) {
      ret <- model.tree$tree[node, "kernelInstanceID"]
    } else {
      ret <- paste("(",
                   node.to.string(model.tree, ld, show.node.ids),
                   op,
                   node.to.string(model.tree, rd, show.node.ids),
                   ")",
                   sep=" "
      )
    }
  }
  if (show.node.ids) {
    ret <- paste("[", node, ": ", ret, " :", node, "]",
                 sep=""
                 )
  }
  return(ret)
}

#' Convert a Model Tree to a String
#'
#' Return or print a string representation of the combination of kernels represented by a given model tree
#'
#' @param model.tree a modelTree object
#' @param show.node.ids boolean indicating whether the string should be annotated with the associated node numbers.
#'
#' @return a string representing a modelTree object
#' @export
#'
#' @examples
#' mt <- create.model.tree.builtin()
#' mt <- insert.kernel.instance(mt, 1, "SE", NULL)
#' mt <- insert.kernel.instance(mt, 1, "SE", "+")
#'
#' as.character(mt) #"( SE.1 + SE.2 )"
#'
#' as.character(mt, show.node.ids=TRUE) #"[3: ( [1: SE.1 :1] + [2: SE.2 :2] ) :3]"
#'
#' print(mt) #"( SE.1 + SE.2 )"
#'
#' print(mt, show.node.ids=TRUE) #"[3: ( [1: SE.1 :1] + [2: SE.2 :2] ) :3]"
as.character.ModelTree <- function(model.tree, show.node.ids=FALSE) {
  root.node <- find.root.node(model.tree)
  return(node.to.string(model.tree, root.node, show.node.ids))
}

#' @rdname as.character.ModelTree
#' @export
print.ModelTree <- function(model.tree, show.node.ids=FALSE) {
  cat(as.character.ModelTree(model.tree, show.node.ids))
  cat("\n")
}

#' Print a Summary of a Model Tree
#'
#' @param model.tree ModelTree object to summarise
#'
#' @export
#'
summary.ModelTree <- function(model.tree) {
  cat("Operation functions:\n")
  for (op in names(model.tree$operation.functions)) {
    cat(paste('- "', op, '"\n', sep=""))
  }

  cat("\nKernel Classes:\n")
  for (kernel in names(model.tree$kernel.objects)) {
    cat(paste('- ', kernel, '\n', sep=""))
    print(model.tree$kernel.objects[[kernel]])
  }

  cat("\nKernel Instances:\n")
  for (i in seq_along(model.tree$kernel.instances[, 1])) {
    kinst_ID <- model.tree$kernel.instances$kernelInstanceID[i]
    kinst_class <- model.tree$kernel.instances$kernelClassName[i]
    cat(paste('- ', kinst_ID, '\n', sep=""))
    cat(paste('  - Class:', kinst_class, '\n', sep=""))
    cat(paste('  - Hyperparameters:\n', sep=""))
    hp_vec <- model.tree$kernel.inst.hyper.param.indices[[kinst_ID]]
    for (hp in seq_along(hp_vec)) {
      if (is.null(model.tree$all.hyper.params[hp_vec[hp]])) {
        hp_val <- "NULL"
      } else {
        hp_val <- model.tree$all.hyper.params[hp_vec[hp]]
      }
      cat(paste('    - ', model.tree$kernel.objects[[kinst_class]][["hyperparam_names"]][hp],
                ": ", hp_val, '\n', sep=""))
    }
  }

  cat("\nModel Tree Formula:\n")
  cat(paste("  ", as.character(model.tree), "\n", sep=""))

  cat("\nModel Tree Formula (with node number annotations):\n")
  cat(paste("  ", as.character(model.tree, show.node.ids=T), "\n", sep=""))
}

#' Add a Kernel Class to a Model Tree
#'
#' @param orig.model.tree a ModelTree object
#' @param kernel.class.name the name of the kernel class (arbitrary, must not
#' match the class name of any kernel already in the model tree).
#' @param kernel either an R function defining a kernel, or a string indicating a built-in kernel type
#' @param hyper.param.names a vector containing the names of the kernel hyperparameters
#' @param kernel.additional.params a list containing any additional kernel parameters
#' @param kernel.grad either an R function returning the gradient of \code{kernel}, or NULL. If \code{kernel}
#' is a string specifying a built-in kernel, \code{kernel.grad} should be NULL. If \code{kernel} is an
#' R function and \code{kernel.grad} is NULL, \code{\link[numDeriv]{grad}} will be used to derive a numerical
#' gradient.
#'
#' @return a ModelTree object
#' @export
#'
#' @examples
#' mt <- create.model.tree()
#' mt <- add.kernel(mt, "SE", "squaredExponential", c("l"))
add.kernel <- function(orig.model.tree,
                       kernel.class.name,
                       kernel) {
  model.tree <- orig.model.tree

  if (kernel.class.name %in% names(model.tree$kernel.objects)) {
    stop("Kernel with that name is already included")
  }

  if (grepl(pattern="#", x=kernel.class.name)) {
    stop("Kernel class name cannot contain '#'")
  }

  model.tree$kernel.objects[[kernel.class.name]] <- kernel


  return(model.tree)
}

#' Insert Kernel Instance into Model Tree
#'
#' @param orig.model.tree a ModelTree object
#' @param node the node at which to insert the kernel instance
#' @param kernel.class.name a string indicating the class of kernel to insert
#' @param operation.id a string indicating how to insert the kernel at \code{node}.
#' Must be "*" or "+", or NULL if the model tree is currently empty.
#' @param hyper.params the hyper parameters of the inserted kernel. Can be NULL to
#' fit hyper parameters at a later time.
#'
#' @return a ModelTree object
#' @export
#'
#' @examples
#' mt <- create.model.tree.builtin()
#'
#' mt <- insert.kernel.instance(mt, 1, "SE", NULL)
#' print(mt) # "SE.1"
#'
#' mt <- insert.kernel.instance(mt, 1, "SE", "*")
#' print(mt) # "( SE.1 * SE.2 )"
#'
#' mt <- insert.kernel.instance(mt, 3, "RQ", "+")
#' print(mt) # "( ( SE.1 * SE.2 ) + RQ.1 )"
#'
insert.kernel.instance <- function(orig.model.tree,
                                   node,
                                   kernel.class.name,
                                   operation.id,
                                   hyper.params=NULL) {
  model.tree <- orig.model.tree

  if (!kernel.class.name %in% names(model.tree$kernel.objects)) {
    stop("Unknown kernel specified.")
  }
  if (node > max(nrow(model.tree$tree), 1) | node <= 0 | floor(node) != node) {
    stop("Invalid node specified.")
  }

  # Check there's no operation if it's the first node
  if(nrow(model.tree$tree)==0 & node==1) {
    if (!is.null(operation.id)) {
      stop("Tree is empty - no operation should be specified for first node inserted.")
    }
  } else {
    # If it's not the first node we need a valid operation.
    if (is.null(operation.id)) {
      stop("No operation specified.")
    } else if (!operation.id %in% names(model.tree$operation.functions)) {
      stop("Invalid operation specified")
    }
  }

  # Check the hyperparameters look right
  hyperparam_names <- model.tree$kernel.objects[[kernel.class.name]]$hyperparam_names
  if (length(hyperparam_names) > 0) {
    if (!is.null(hyper.params)) {
      # If we expect some and have been given some, do the two match?
      if (!all(names(hyper.params) %in% hyperparam_names)) {
        stop("Supplied hyperparameters do not match hyperparameters expected for this kernel.")
      }
    }
  } else {
    # We're not expecting any hyperparams.
    if (!is.null(hyper.params)) {
      stop("Hyperparameters supplied when none are expected for this kernel.")
    }
  }


  inst.id <- generate.id(model.tree$kernel.instances[1,], kernel.class.name)

  # Add to the list of instances
  model.tree$kernel.instances[inst.id, ] <- c(kernelInstanceID=inst.id,
                                              kernelClassName=kernel.class.name
  )


  # Append or create hyper parameters. NAs indicate random start.
  # Are we expecting any hyperparameters?
  if (length(hyperparam_names) > 0) {
    # Have any been supplied?
    if (!is.null(hyper.params)) {
      # Re-order if necessary to match our stored list.
      hyper.params <- hyper.params[hyperparam_names]
    } else {
      # Make NA hyperparams to indicate random start
      hyper.params <- rep(NA, length(hyperparam_names))
    }
    names(hyper.params) <- paste(inst.id, hyperparam_names, sep="#")

    # Add them to the list of hyper.params
    start.index <- length(model.tree$all.hyper.params) + 1
    model.tree$all.hyper.params <- c(model.tree$all.hyper.params, hyper.params)
    end.index <- length(model.tree$all.hyper.params)

    model.tree$kernel.inst.hyper.param.indices[[inst.id]] <- start.index:end.index

  } else {
    # We're not expecting any hyperparams. Create an empty hyperparam indices list.
    model.tree$kernel.inst.hyper.param.indices[[inst.id]] <- c()
  }

  kernel.inst.node.entry <- c(leftDaughter=NA,
                              rightDaughter=NA,
                              kernelInstanceID=inst.id,
                              operationID=NA
  )

  # Find the root node of this model tree. environment() returns the containing model tree env.
  root.node <- find.root.node(model.tree)

  # Are we inserting the first node?
  if (nrow(model.tree$tree)==0 & node==1) {
    model.tree$tree[1, ] <- kernel.inst.node.entry
  } else {
    new.inst.node <- nrow(model.tree$tree) + 1
    new.op.node <- nrow(model.tree$tree) + 2


    new.op.node.entry <- c(leftDaughter=node,
                           rightDaughter=new.inst.node,
                           kernelInstanceID=NA,
                           operationID=operation.id
    )

    # If node is the root node we're done. Otherwise we need to fix the parent
    # to point at our new operation node.
    if (node!=root.node) {
      # Find the parent to this node so we can have it point at our new
      # operation node rather than the current node.
      parent.node <- find.parent.node(model.tree, node)

      node.is.left.daughter <- model.tree$tree[parent.node, "leftDaughter"] == node

      # point the parent at our new op node
      if (node.is.left.daughter) {
        model.tree$tree[parent.node, "leftDaughter"] <- new.op.node
      } else {
        # If it's not left then it must be right
        model.tree$tree[parent.node, "rightDaughter"] <- new.op.node
      }
    }

    # Insert our new nodes
    model.tree$tree[new.inst.node, ] <- kernel.inst.node.entry
    model.tree$tree[new.op.node, ] <- new.op.node.entry

  }
  return(model.tree)
}

#' Delete a Node in a Model Tree
#'
#' @param model.tree a ModelTree object
#' @param node the node to delete (must not be the root node)
#'
#' @return a copy of \code{model.tree} with \code{node} deleted.
#' @export
delete.node <- function(model.tree, node) {
  root.node <- find.root.node(model.tree)
  if (node == root.node) {
    stop("Can't delete the root node.")
  }

  new.mt <- model.tree

  # Repoint the grandparent node at this node's sibling
  parent.node <- find.parent.node(new.mt, node)

  if (new.mt$tree[parent.node, "leftDaughter"] == node) {
    sibling.node <- new.mt$tree[parent.node, "rightDaughter"]
    new.mt$tree[parent.node, "rightDaughter"] <- NA
  } else {
    sibling.node <- new.mt$tree[parent.node, "leftDaughter"]
    new.mt$tree[parent.node, "leftDaughter"] <- NA
  }

  if (!parent.node == root.node) {
    grandparent.node <- find.parent.node(new.mt, parent.node)

    if (new.mt$tree[grandparent.node, "leftDaughter"] == parent.node) {
      new.mt$tree[grandparent.node, "leftDaughter"] <- sibling.node
    } else {
      new.mt$tree[grandparent.node, "rightDaughter"] <- sibling.node
    }
  }

  # Now delete the orphaned nodes (by overwriting with NAs then deleting the rows)
  deleted.nodes <- numeric(0)
  node.queue <- c(parent.node)
  while (length(node.queue) > 0) {
    current.node <- node.queue[1]
    node.queue <- node.queue[-1]
    node.row <- new.mt$tree[current.node, ]
    if (!is.na(node.row$operationID)) {
      node.queue <- as.numeric(c(node.queue,
                      new.mt$tree[current.node, "leftDaughter"],
                      new.mt$tree[current.node, "rightDaughter"]))
      node.queue <- node.queue[!is.na(node.queue)]
    }
    new.mt$tree[current.node, ] <- NA

    deleted.nodes <- c(deleted.nodes, current.node)
  }

  node.indices <- seq(nrow(new.mt$tree))[-deleted.nodes]

  new.mt$tree <- new.mt$tree[-deleted.nodes, ]

  rownames(new.mt$tree) <- 1:nrow(new.mt$tree)

  # We need to update the node references in the tree
  for (i in 1:nrow(new.mt$tree)) {
    for (child in c("leftDaughter", "rightDaughter")) {
      if (!is.na(new.mt$tree[i, child])) {
        new.mt$tree[i, child] <- which(node.indices == new.mt$tree[i, child])
      }
    }
  }

  remaining.inst.ids <- unique(new.mt$tree[, "kernelInstanceID"])
  deleted.hyper.param.indices <- c()

  for (i in seq_along(new.mt$kernel.instances[, "kernelInstanceID"])) {
    inst.id <- new.mt$kernel.instances[i, "kernelInstanceID"]
    if (!(inst.id %in% remaining.inst.ids)) {
      new.mt$kernel.instances <- new.mt$kernel.instances[-i, ]

      deleted.hyper.param.indices <- c(deleted.hyper.param.indices,
                                       new.mt$kernel.inst.hyper.param.indices[[inst.id]])

      temp.index <- which(names(new.mt$kernel.inst.hyper.param.indices) == inst.id)
      new.mt$kernel.inst.hyper.param.indices <- new.mt$kernel.inst.hyper.param.indices[-temp.index]
    }
  }

  # If we deleted some hyper parameters we need to translate the old indices to
  # the new indices.
  if (length(deleted.hyper.param.indices) > 0) {
    hyper.param.indices <- seq_along(new.mt$all.hyper.params)[-deleted.hyper.param.indices]
    new.mt$all.hyper.params <- new.mt$all.hyper.params[-deleted.hyper.param.indices]
    for (i in seq_along(new.mt$kernel.inst.hyper.param.indices)) {
      for (j  in seq_along(new.mt$kernel.inst.hyper.param.indices[[i]])) {
        new.mt$kernel.inst.hyper.param.indices[[i]][j] <-
          which(hyper.param.indices == new.mt$kernel.inst.hyper.param.indices[[i]][j])
      }
    }
  }

  return(new.mt)
}

generate.next.models <- function(model.tree) {

  # Generate a list of models
  models <- list()
  i <- 1
  if (nrow(model.tree$tree) == 0) {
    for (kernel in names(model.tree$kernel.objects)) {
      models[[i]] <- insert.kernel.instance(orig.model.tree=model.tree,
                                            node=1,
                                            kernel.class.name=kernel,
                                            operation.id=NULL)
      i <- i + 1
    }
  } else {
    for (operation in names(model.tree$operation.functions)) {
      for (kernel in names(model.tree$kernel.objects)) {
        for (node in 1:nrow(model.tree$tree)) {
          models[[i]] <- insert.kernel.instance(orig.model.tree=model.tree,
                                                node=node,
                                                kernel.class.name=kernel,
                                                operation.id=operation)
          i <- i + 1
        }
      }
    }
    root.node <- find.root.node(model.tree)
    for (node in 1:nrow(model.tree$tree)) {
      if (node != root.node) {
        models[[i]] <- delete.node(model.tree, node)
        i <- i + 1
      }
    }
  }

  # Reduce to canonical form and dedupe

  model.strings <- character(length(models))
  for (i in seq_along(models)) {
    models[[i]] <- reduce.to.canonical.tree(models[[i]])
    model.strings[i] <- as.character(models[[i]])
  }

  duplicates <- which(duplicated(model.strings))

  if (length(duplicates) > 0) {
    models <- models[-duplicates]
  }

  return(models)
}

#' Create a Model Tree from a list of kernel objects
#'
#' Returns a model tree instance populated with a basis set of kernels.
#'
#' @param kernel_list a named list of kernel objects
#'
#' @return A modelTree object
#' @export
create.model.tree.from.list <- function(kernel_list) {

  mt <- create.model.tree()
  for (kernel_name in names(kernel_list)) {
    kernel <- kernel_list[[kernel_name]]
    mt <- add.kernel(mt,
                     kernel_name,
                     kernel
    )
  }

  return(mt)
}

#' Create a Model Tree Containing All Built-in Kernels
#'
#' Returns a model tree instance populated with a basis set of all of the built in kernels.
#'
#' @return A modelTree object
#' @export
#'
#' @examples
#' mt <- create.model.tree.builtin()
create.model.tree.builtin <- function() {
  builtin.kernels <- list_built_in_kernels()

  mt <- create.model.tree.from.list(builtin.kernels)

  return(mt)
}

#' Create a Model Tree Containing All Built-in Kernels (scaled)
#'
#' Returns a model tree instance populated with a basis set of all of the built in kernels.
#'
#' @return A modelTree object
#' @export
#'
#' @examples
#' mt <- create.model.tree.builtin.scaled()
create.model.tree.builtin.scaled <- function() {
  builtin.kernels <- list_built_in_kernels_scaled()

  mt <- create.model.tree.from.list(builtin.kernels)

  return(mt)
}


check.node.is.in.canonical.form <- function(model.tree, node) {
  is.canonical <- TRUE

  node.row <- model.tree$tree[node, ]

  # If this is a terminal (non-operation) node then it is canonical. If not...
  if (!is.na(node.row$operationID)) {
    # check left daughter operation does not match this node's
    if (!identical(model.tree$tree[node.row$leftDaughter, "operationID"], node.row$operationID)) {
      is.canonical <- FALSE
    } else {
      # check the left daughter of this node is smaller than the left daughter of it's right daughter,
      # if they are the same type of operation
      right.daughter.node.row <- model.tree$tree[node.row$rightDaughter, ]
      if (identical(right.daughter.node.row$operationID, node.row$operationID)) {
        right.daughters.left.daughter.string <- node.to.string(right.daughter.node.row$leftDaughter)
        left.daughter.string <- node.to.string(node.row$leftDaughter)
        if (right.daughters.left.daughter.string < left.daughter.string) {
          is.canonical <- FALSE
        }
      }
    }
  }
  return(is.canonical)
}

block.of.ops.is.stretched <- function(model.tree, node, operationID) {
  node.row <- model.tree$tree[node, ]

  if (is.na(operationID)) {
    return(TRUE)
  }

  if (identical(operationID, node.row$operationID)){
    left.daughter.op <- model.tree$tree[node.row$leftDaughter, "operationID"]
    if (identical(operationID, left.daughter.op)) {
      return(FALSE)
    } else {
      return(block.of.ops.is.stretched(model.tree, node.row$rightDaughter, operationID))
    }
  } else {
    return(TRUE)
  }
}

find.block.terminal.node <- function(model.tree, node, operationID) {
  tree <- model.tree$tree
  node.row <- tree[node,]
  if (!identical(node.row$operationID, operationID)) {
    stop("block root node is not of the specified operation")
  }
  current.node <- node
  next.node <- tree[current.node, "rightDaughter"]
  while (identical(tree[next.node, "operationID"], operationID)) {
    current.node <- next.node
    next.node <- tree[next.node, "rightDaughter"]
  }
  return(current.node)
}

stretch.block <- function(model.tree, node, operationID) {
  node.row <- model.tree$tree[node,]
  if (!identical(node.row$operationID, operationID)) {
    stop("block root node is not of the specified operation")
  }
  current.node <- node
  terminal.node <- find.block.terminal.node(model.tree, current.node, operationID)
  while (identical(model.tree$tree[current.node, "operationID"], operationID)) {
    left.daughter <- model.tree$tree[current.node, "leftDaughter"]
    if (identical(model.tree$tree[left.daughter, "operationID"], operationID)) {
      # We need a non-matching op as left daughter, so we'll swap the left
      # daughter with the right daughter of the terminal node
      terminal.node.right.daughter <- model.tree$tree[terminal.node, "rightDaughter"]
      model.tree$tree[current.node, "leftDaughter"] <- terminal.node.right.daughter
      model.tree$tree[terminal.node, "rightDaughter"] <- left.daughter

      # Now we'll have a new terminal node, so let's go find it
      terminal.node <- find.block.terminal.node(model.tree, current.node, operationID)
    }

    # We're done with this node, so let's move on.
    current.node <- model.tree$tree[current.node, "rightDaughter"]
  }
  return(model.tree)
}

stretch.tree <- function(model.tree) {
  root.node <- find.root.node(model.tree)
  node.queue <- c(root.node)
  while (length(node.queue) > 0) {
    #Pop off of the front of our wildly inefficient "queue"
    current.node <- node.queue[1]
    node.queue <- node.queue[-1]

    current.op <- model.tree$tree[current.node, "operationID"]

    if (!block.of.ops.is.stretched(model.tree, current.node, current.op)) {
      model.tree <- stretch.block(model.tree, current.node, current.op)
    }

    if (!is.na(current.op)) {
      node.queue <- c(node.queue, as.numeric(model.tree$tree[current.node, c("leftDaughter", "rightDaughter")]))
    }
  }
  return(model.tree)
}

order.and.return.string <- function(model.tree, node) {
  # We assume the tree has been stretched i.e. commutable operations proceed
  # in blocks down the right daughters
  node.op <- model.tree$tree[node, "operationID"]
  if (!is.na(node.op)) {
    left.daughter <- model.tree$tree[node, "leftDaughter"]
    right.daughter <- model.tree$tree[node, "rightDaughter"]
    right.daughter.op <- model.tree$tree[right.daughter, "operationID"]
    if (!identical(right.daughter.op, node.op)) {
      # Terminal string of a block - get strings and re-order if necessary
      left.temp <- order.and.return.string(model.tree, left.daughter)
      left.daughter.string <- left.temp$string
      model.tree <- left.temp$model.tree

      right.temp <- order.and.return.string(model.tree, right.daughter)
      right.daughter.string <- right.temp$string
      model.tree <- right.temp$model.tree
      if (left.daughter.string > right.daughter.string) {
        model.tree$tree[node, "leftDaughter"] <- right.daughter
        model.tree$tree[node, "rightDaughter"] <- left.daughter
      }
    } else {
      # We're in a block of ops - need to look at the ops around us.
      # If we're in a block then "node" must be the root node of the
      # block, since we never call this function on right daughters in
      # the block.

      # Get the block's left daughter strings
      current.block.node <- node
      block.nodes <- c()
      block.left.daughters <- c()
      block.left.daughter.strings <- c()
      while (identical(model.tree$tree[current.block.node, "operationID"], node.op)) {
        block.nodes <- c(block.nodes, current.block.node)
        block.left.daughter <- model.tree$tree[current.block.node, "leftDaughter"]
        block.left.daughters <- c(block.left.daughters, block.left.daughter)
        block.left.daughter.temp <- order.and.return.string(model.tree,
                                                            block.left.daughter)
        model.tree <- block.left.daughter.temp$model.tree
        block.left.daughter.strings <- c(block.left.daughter.strings,
                                         block.left.daughter.temp$string)
        current.block.node <- model.tree$tree[current.block.node, "rightDaughter"]
      }

      # current.block.node now contains the right daughter of the terminal node -
      # we need to include this in the ordering
      terminal.right.daughter <- current.block.node
      block.left.daughters <- c(block.left.daughters, terminal.right.daughter)
      block.left.daughter.temp <- order.and.return.string(model.tree,
                                                          terminal.right.daughter)
      model.tree <- block.left.daughter.temp$model.tree
      block.left.daughter.strings <- c(block.left.daughter.strings,
                                       block.left.daughter.temp$string)

      correct.daughter.order <- order(block.left.daughter.strings)
      ordered.daughters <- block.left.daughters[correct.daughter.order]
      for (block.node.index in seq_along(block.nodes)) {
        block.node <- block.nodes[block.node.index]
        new.left.daughter <- ordered.daughters[block.node.index]
        model.tree$tree[block.node, "leftDaughter"] <- new.left.daughter
      }

      # Now we update the right daughter of the terminal node
      block.node <- tail(block.nodes, 1)
      new.right.daughter <- tail(ordered.daughters, 1)
      model.tree$tree[block.node, "rightDaughter"] <- new.right.daughter
    }
  }
  return(list(string=node.to.string(model.tree, node, show.node.ids=FALSE),
              model.tree=model.tree))
}

order.tree <- function(model.tree) {
  root.node <- find.root.node(model.tree)
  model.temp <- order.and.return.string(model.tree, root.node)
  model.tree <- model.temp$model.tree
  return(model.tree)
}

reduce.to.canonical.tree <- function(model.tree) {
  # First we stretch the tree, then we order it.
  model.tree <- stretch.tree(model.tree)
  model.tree <- order.tree(model.tree)
  return(model.tree)
}
