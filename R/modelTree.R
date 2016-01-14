clone.env <- function(orig.env) {
  copy.env <- new.env()
  for(n in ls(orig.env, all.names=TRUE)) {
    assign(n, get(n, orig.env), copy.env)
  }
  attributes(copy.env) <- attributes(orig.env)
  return(copy.env)
}

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

create.model.tree <- function() {
  model.tree <- new.env()
  class(model.tree) <- "ModelTree"

  model.tree$operation.functions = list("+"=kern.sum,
                                        "*"=kern.prod
  )

  model.tree$operation.grad.functions = list("+"=sum.grad,
                                             "*"=prod.grad
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

  # This holds the functions (or string IDs for built-in kernels) of kernels
  # available to the model. Indexed by kernelClassName.
  model.tree$kernel.class.functions <- list()
  model.tree$kernel.class.grad.functions <- list()
  model.tree$kernel.class.additional.params <- list()

  # This holds the hyper parameter names for each kernel. Indexed by kernelClassName
  model.tree$kernel.class.hyper.param.names <- list()

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

as.character.ModelTree <- function(model.tree, show.node.ids=FALSE) {
  root.node <- find.root.node(model.tree)
  return(node.to.string(model.tree, root.node, show.node.ids))
}

print.ModelTree <- function(model.tree, show.node.ids=FALSE) {
  print(as.character.ModelTree(model.tree, show.node.ids))
}

add.kernel <- function(orig.model.tree,
                       kernel.class.name,
                       kernel,
                       hyper.param.names,
                       kernel.additional.params=list(),
                       kernel.grad=NULL) {
  model.tree <- clone.env(orig.model.tree)

  if (kernel.class.name %in% names(model.tree$kernel.class.functions)) {
    stop("Kernel with that name is already included")
  }
  model.tree$kernel.class.functions[[kernel.class.name]] <- kernel
  model.tree$kernel.class.additional.params[[kernel.class.name]] <- kernel.additional.params

  if (is.null(kernel.grad)) {
    if (is.character(kernel)) {
      # Built-in kernel
      model.tree$kernel.class.grad.functions[[kernel.class.name]] <- kernel
    } else {
      # No grad supplied. Make a grad.
      model.tree$kernel.class.grad.functions[[kernel.class.name]] <- create.numeric.grad(kernel)
    }
  } else {
    # Grad supplied.
    model.tree$kernel.class.grad.functions[[kernel.class.name]] <- kernel.grad
  }

  model.tree$kernel.class.hyper.param.names[[kernel.class.name]] <- hyper.param.names

  return(model.tree)
}

insert.kernel.instance <- function(orig.model.tree,
                                   node,
                                   kernel.class.name,
                                   operation.id,
                                   hyper.params=NULL) {
  model.tree <- clone.env(orig.model.tree)

  if (!kernel.class.name %in% names(model.tree$kernel.class.functions)) {
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
  if (length(model.tree$kernel.class.hyper.param.names[[kernel.class.name]]) > 0) {
    if (!is.null(hyper.params)) {
      # If we expect some and have been given some, do the two match?
      if (!all(names(hyper.params) %in% model.tree$kernel.class.hyper.param.names[[kernel.class.name]])) {
        stop("Supplied hyperparameters do not match hyperparameters expected for this kernel.")
      }
    }
  } else {
    # We're not expecting any hyperparams.
    if (!is.null(hyper.params)) {
      stop("Hyperparameters supplied when none are expected for this kernel.")
    }
  }


  inst.id <- generate.id(model.tree$tree$kernelInstanceID, kernel.class.name)

  # Add to the list of instances
  model.tree$kernel.instances[inst.id, ] <- c(kernelInstanceID=inst.id,
                                              kernelClassName=kernel.class.name
  )


  # Append or create hyper parameters. NAs indicate random start.
  # Are we expecting any hyperparameters?
  if (length(model.tree$kernel.class.hyper.param.names[[kernel.class.name]]) > 0) {
    # Have any been supplied?
    if (!is.null(hyper.params)) {
      # Re-order if necessary to match our stored list.
      hyper.params <- hyper.params[model.tree$kernel.class.hyper.param.names[[kernel.class.name]]]
    } else {
      # Make NA hyperparams to indicate random start
      hyper.params <- rep(NA, length(model.tree$kernel.class.hyper.param.names[[kernel.class.name]]))
      names(hyper.params) <- model.tree$kernel.class.hyper.param.names[[kernel.class.name]]
    }

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

generate.next.models <- function(model.tree) {
  models <- list()
  i <- 1
  if (nrow(model.tree$tree) == 0) {
    for (kernel in names(model.tree$kernel.class.functions)) {
      models[[i]] <- insert.kernel.instance(orig.model.tree=model.tree,
                                          node=1,
                                          kernel.class.name=kernel,
                                          operation.id=NULL)
      i <- i + 1
    }
  } else {
    for (operation in names(model.tree$operation.functions)) {
      for (kernel in names(model.tree$kernel.class.functions)) {
        for (node in 1:nrow(model.tree$tree)) {
          models[[i]] <- insert.kernel.instance(orig.model.tree=model.tree,
                                                node=node,
                                                kernel.class.name=kernel,
                                                operation.id=operation)
          i <- i + 1
        }
      }
    }
  }
  return(models)
}

create.model.tree.builtin <- function() {
  mt <- create.model.tree()
  for (kernel in builtin.kernels) {
    mt <- add.kernel(mt,
                     kernel,
                     kernel,
                     builtin.kernel.hyperparam.names[[kernel]],
                     builtin.kernel.additional.hyperparams[[kernel]]
                     )
  }
  return(mt)
}

reduce.to.canonical.tree <- function(model.tree) {
  # kernel instances are ordered lexicographically by ID.
  # First, any products are multiplied out so (a + b) * c becomes a * c + b * c,
  # a * (b + c) become a * b + a * c, and
}
