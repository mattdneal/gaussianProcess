
ahp.cache.name <- "all.hyper.params.cache"
icm.cache.name <- "inst.cov.mats.cache"
icmg.cache.name <- "inst.cov.mats.grad.cache"


#' Clone an Environment
#'
#' @param orig.env the environment to clone
#'
#' @return a new environment whose contents match \code{orig.env}
#' @export
#'
#' @seealso \code{\link{clone.ModelTree}}
#'
#' @examples
#' mt <- create.model.tree.builtin()
#' mt2 <- clone.env(mt)
clone.env <- function(orig.env) {
  copy.env <- new.env()
  for(n in ls(orig.env, all.names=TRUE)) {
    assign(n, get(n, orig.env), copy.env)
  }
  attributes(copy.env) <- attributes(orig.env)
  return(copy.env)
}

#' Clone a Model Tree
#'
#' @param orig.model.tree the ModelTree to clone
#' @param invalidate.cache boolean indicating whether the cache should be cleared.
#' If you're not sure whether you need to clear the cache or not, you probably do.
#'
#' @return a new Model Tree whose contents match \code{orig.model.tree}
#' @export
#'
#' @examples
#' mt <- create.model.tree.builtin()
#' mt2 <- clone.ModelTree(mt)
clone.ModelTree <- function(orig.model.tree, invalidate.cache=TRUE) {
  model.tree <- clone.env(orig.model.tree)
  if (invalidate.cache) {
    invalidate.cache(model.tree)
  }
  return(model.tree)
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

#' Create an Empty Model Tree
#'
#' Creates an instance of a modelTree object with no associated kernels.
#'
#' N.B. ModelTree objects are implemented using environments, not lists,
#' and as such need to be explicitly cloned using \code{\link{clone.env}}.
#'
#' @return a modelTree object
#' @export
#'
#' @examples
#' mt <- create.model.tree()
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

  # Calling invalidate.cache will set up an empty cache for us
  invalidate.cache(model.tree)

  return(model.tree)
}

# Because model trees are environments we can update the cache and it persists
# outside of the function.
invalidate.cache <- function(model.tree) {
  model.tree$cache <- list()
  model.tree$cache[[icm.cache.name]] <- list()
  model.tree$cache[[icmg.cache.name]] <- list()
  return(NULL)
}

# Because model trees are environments we can update the cache and it persists
# outside of the function.
hyper.param.cache.check <- function(model.tree) {
  if (ahp.cache.name %in% names(model.tree$cache)) {
    ahp.cache <- model.tree$cache[[ahp.cache.name]]
    ahp <- model.tree$all.hyper.params

    for (kinst.ID in model.tree$kernel.instances$kernelInstanceID) {
      hp.indices <- model.tree$kernel.inst.hyper.param.indices[[kinst.ID]]

      if (any(ahp[hp.indices] != ahp.cache[hp.indices])) {
        model.tree$cache[[icm.cache.name]][[kinst.ID]] <- NULL
        model.tree$cache[[icmg.cache.name]][[kinst.ID]] <- NULL
      }
    }
  }

  model.tree$cache[[ahp.cache.name]] <- model.tree$all.hyper.params
  return(NULL)
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
  for (kclass in names(model.tree$kernel.class.functions)) {
    cat(paste('- ', kclass, '\n', sep=""))
    if (is.character(model.tree$kernel.class.functions[[kclass]])) {
      cat(paste('  - Built-in kernel: ', model.tree$kernel.class.functions[[kclass]], '\n', sep=""))
    }
    if (length(model.tree$kernel.class.additional.params[[kclass]]) > 0) {
      cat(paste('  - Additional Parameters:\n', sep=""))
      for (addparam in names(model.tree$kernel.class.additional.params[[kclass]])) {
        cat(paste('    - ', addparam, '\n', sep=""))
      }
    }
  }

  cat("\nKernel Instances:\n")
  for (i in 1:nrow(model.tree$kernel.instances)) {
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
      cat(paste('    - ', model.tree$kernel.class.hyper.param.names[[kinst_class]][hp],
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
                       kernel,
                       hyper.param.names,
                       kernel.additional.params=list(),
                       kernel.grad=NULL) {
  model.tree <- clone.ModelTree(orig.model.tree)

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
  model.tree <- clone.ModelTree(orig.model.tree)

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

#' Add Polynomial Kernel Classes to a Model Tree
#'
#' @param mt a ModelTree object
#' @param degrees a vector containing the degrees to add
#' @param homogeneous boolean indicating whether to add homogeneous polynomial kernels
#' @param poly boolean indicating whether to add non-homogeneous polynomial kernels
#' @param generalised boolean indicating whether to add generealised polynomial kernels
#'
#' @return a ModelTree object
#' @export
#'
#' @examples
#' mt <- create.model.tree()
#' mt <- add.polynomial.kernels(mt, degrees=c(2,3))
add.polynomial.kernels <- function(mt, degrees=c(1,2),
                                   homogeneous=TRUE, poly=TRUE, generalised=TRUE) {
  for (degree in degrees) {
    if (homogeneous) {
      mt <- add.kernel(mt,
                       paste("homogPoly", degree, sep="_"),
                       "homogeneousPolynomial",
                       numeric(0),
                       list(degree=degree)
      )
    }
    if (poly) {
      mt <- add.kernel(mt,
                       paste("poly", degree, sep="_"),
                       "polynomial",
                       c("c"),
                       list(degree=degree)
      )
    }
    if (generalised) {
      mt <- add.kernel(mt,
                       paste("genPoly", degree, sep="_"),
                       "generalisedPolynomial",
                       c("l", "c"),
                       list(degree=degree)
      )
    }
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
  builtin.kernels <- c(SE="squaredExponential",
                       RQ="rationalQuadratic",
                       periodic="periodic",
                       constant="constant",
                       NN="neuralNetwork"

  )

  # If there are no hyperparams, use numeric(0) rather than c() so the c++
  # code works - it expects a vector not NULL.
  builtin.kernel.hyperparam.names <- list(SE=c("l"),
                                          RQ=c("l", "alpha"),
                                          periodic=c("l", "p"),
                                          constant=c("sigma_0"),
                                          NN=c("sigma_0", "sigma")
                                          #"oneDLinear"=c("intercept","sigma_1")
                                          #"changepoint"=c("changepoint", "transitionRate")
  )

  builtin.kernel.additional.hyperparams <- list(SE=list(),
                                                RQ=list(),
                                                periodic=list(),
                                                constant=list(),
                                                NN=list()
                                                #"oneDLinear"=list(dimensionIndex=1)
                                                #"changepoint"=list(dimensionIndex=1)
  )

  mt <- create.model.tree()
  for (kernel_index in seq_along(builtin.kernels)) {
    kernel_name <- names(builtin.kernels)[kernel_index]
    mt <- add.kernel(mt,
                     kernel_name,
                     builtin.kernels[kernel_name],
                     builtin.kernel.hyperparam.names[[kernel_name]],
                     builtin.kernel.additional.hyperparams[[kernel_name]]
                     )
  }

  mt <- add.polynomial.kernels(mt, degrees=c(1,2), homogeneous=FALSE, poly=FALSE)

  return(mt)
}

reduce.to.canonical.tree <- function(model.tree) {
  # kernel instances are ordered lexicographically by ID.
  # First, any products are multiplied out so (a + b) * c becomes a * c + b * c,
  # a * (b + c) become a * b + a * c, and
}
