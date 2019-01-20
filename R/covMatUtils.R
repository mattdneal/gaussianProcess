#' Evaluate Covariance Matrix of a Model Tree Node
#'
#'  This function recursively evaluates a model tree at the given node,
#'  returning the covariance matrix for that node
#'
#' @param model.tree an object of class modelTree
#' @param inst.cov.mats a named list of the covariance matrices of the kernel instances in the model tree
#' @param node the node to evaluate
#'
#' @return A covariance matrix
#'
eval.cov.mat.node <- function(model.tree, inst.cov.mats, node) {
  node.entry <- model.tree$tree[node, ]
  if (!is.na(node.entry["operationID"])) {
    opID <- as.character(node.entry["operationID"])
    operation <- model.tree$operation.functions[[opID]]

    left.daughter <- as.numeric(node.entry["leftDaughter"])
    right.daughter <- as.numeric(node.entry["rightDaughter"])

    left.daughter.mat <- eval.cov.mat.node(model.tree, inst.cov.mats, left.daughter)
    right.daughter.mat <- eval.cov.mat.node(model.tree, inst.cov.mats, right.daughter)

    return(operation(left.daughter.mat, right.daughter.mat))
  } else if (!is.na(node.entry["kernelInstanceID"])) {
    inst.id <- as.character(node.entry["kernelInstanceID"])
    return(inst.cov.mats[[inst.id]])
  } else {
    stop(paste("operationID and kernelInstanceID both NULL for node", node))
  }
}

#' Evaluate Covariance Matrix Grad of a Model Tree Node
#'
#'  This function recursively evaluates a model tree at the given node, returning
#'  the grad of the covariance matrix at that node
#'
#' @param model.tree An object of class modelTree
#' @param inst.cov.mats a named list of the covariance matrices of the kernel instances in the model tree
#' @param inst.cov.mat.grads The covariance matrix grads of the kernel instances in the model tree
#' @param node The node to evaluate
#'
#' @return A covariance matrix grad (as an array)
eval.cov.mat.grad.node <- function(model.tree, inst.cov.mats, inst.cov.mat.grads, node) {
  node.entry <- model.tree$tree[node, ]
  if (!is.na(node.entry["operationID"])) {
    opID <- as.character(node.entry["operationID"])
    operation.grad <- model.tree$operation.grad.functions[[opID]]

    left.daughter <- as.numeric(node.entry["leftDaughter"])
    right.daughter <- as.numeric(node.entry["rightDaughter"])

    operation.grad.args <- names(formals(operation.grad))

    if ("a" %in% operation.grad.args) {
      a <- eval.cov.mat.node(model.tree, inst.cov.mats, left.daughter)
    } else {
      a <- NULL
    }
    if ("b" %in% operation.grad.args) {
      b <- eval.cov.mat.node(model.tree, inst.cov.mats, right.daughter)
    } else {
      b <- NULL
    }
    if ("a.grad" %in% operation.grad.args) {
      a.grad <- eval.cov.mat.grad.node(model.tree, inst.cov.mats, inst.cov.mat.grads, left.daughter)
    } else {
      a.grad <- NULL
    }
    if ("b.grad" %in% operation.grad.args) {
      b.grad <- eval.cov.mat.grad.node(model.tree, inst.cov.mats, inst.cov.mat.grads, right.daughter)
    } else {
      b.grad <- NULL
    }

    return(operation.grad(a=a, b=b, a.grad=a.grad, b.grad=b.grad))

  } else if (!is.na(node.entry["kernelInstanceID"])) {
    inst.id <- as.character(node.entry["kernelInstanceID"])
    return(inst.cov.mat.grads[[inst.id]])
  } else {
    stop(paste("operationID and kernelInstanceID both NULL for node", node))
  }
}

#' Evaluate Covariance Matrix Hessian of a Model Tree Node
#'
#'  This function recursively evaluates a model tree at the given node, returning
#'  the Hessian of the covariance matrix at that node
#'
#' @param model.tree An object of class modelTree
#' @param inst.cov.mats a named list of the covariance matrices of the kernel instances in the model tree
#' @param inst.cov.mat.grads The covariance matrix grads of the kernel instances in the model tree
#' @param inst.cov.mat.hess The covariance matrix Hessians of the kernel instances in the model tree
#' @param node The node to evaluate
#'
#' @return A covariance matrix grad (as an array)
eval.cov.mat.hess.node <- function(model.tree,
                                   inst.cov.mats,
                                   inst.cov.mat.grads,
                                   inst.cov.mat.hess,
                                   node) {
  node.entry <- model.tree$tree[node, ]
  if (!is.na(node.entry["operationID"])) {
    opID <- as.character(node.entry["operationID"])
    operation.hess <- model.tree$operation.hess.functions[[opID]]

    left.daughter <- as.numeric(node.entry["leftDaughter"])
    right.daughter <- as.numeric(node.entry["rightDaughter"])

    operation.hess.args <- names(formals(operation.hess))

    if ("a" %in% operation.hess.args) {
      a <- eval.cov.mat.node(model.tree, inst.cov.mats, left.daughter)
    } else {
      a <- NULL
    }
    if ("b" %in% operation.hess.args) {
      b <- eval.cov.mat.node(model.tree, inst.cov.mats, right.daughter)
    } else {
      b <- NULL
    }
    if ("a.grad" %in% operation.hess.args) {
      a.grad <- eval.cov.mat.grad.node(model.tree,
                                       inst.cov.mats,
                                       inst.cov.mat.grads,
                                       left.daughter)
    } else {
      a.grad <- NULL
    }
    if ("b.grad" %in% operation.hess.args) {
      b.grad <- eval.cov.mat.grad.node(model.tree,
                                       inst.cov.mats,
                                       inst.cov.mat.grads,
                                       right.daughter)
    } else {
      b.grad <- NULL
    }
    if ("a.hess" %in% operation.hess.args) {
      a.hess <- eval.cov.mat.hess.node(model.tree,
                                       inst.cov.mats,
                                       inst.cov.mat.grads,
                                       inst.cov.mat.hess,
                                       left.daughter)
    } else {
      a.hess <- NULL
    }
    if ("b.hess" %in% operation.hess.args) {
      b.hess <- eval.cov.mat.hess.node(model.tree,
                                       inst.cov.mats,
                                       inst.cov.mat.grads,
                                       inst.cov.mat.hess,
                                       right.daughter)
    } else {
      b.hess <- NULL
    }

    return(operation.hess(a=a, b=b, a.grad=a.grad, b.grad=b.grad, a.hess=a.hess, b.hess=b.hess))

  } else if (!is.na(node.entry["kernelInstanceID"])) {
    inst.id <- as.character(node.entry["kernelInstanceID"])
    return(inst.cov.mat.hess[[inst.id]])
  } else {
    stop(paste("operationID and kernelInstanceID both NULL for node", node))
  }
}


#' Get Kernel Instance Covariance Matrices, Grads or Hessians
#'
#' Calculates the covariance matrices (or associated grads) for the kernel instances appearing in a model tree
#'
#' @param x A numeric matrix
#' @param model.tree A modelTree object
#' @param return.type Whether to return the covariance matrices, the associated grads, or the hessians
#'
#' @importFrom cacheMan cached_call
#'
#' @return A named list containing covariance matrices (or grad arrays), with one entry for each kernel instance.
get.inst.cov.mats <- function(x, model.tree, return.type=c("cov", "grad", "hess"), cache=NULL) {

  return.type <- match.arg(return.type)

  return.grad <- return.type == "grad"

  return.hess <- return.type == "hess"

  inst.cov.mats <- list()

  for (inst.id in rownames(model.tree$kernel.instances)) {

    kernel.class <- model.tree$kernel.instances[inst.id, "kernelClassName"]
    kernel <- model.tree$kernel.objects[[kernel.class]]
    kernel.hyper.param.indices <- model.tree$kernel.inst.hyper.param.indices[[inst.id]]
    kernel.hyper.params <- model.tree$all.hyper.params[kernel.hyper.param.indices]
    model.tree.hp.names <- paste(inst.id, kernel$hyperparam_names, sep="#")
    if (length(kernel$hyperparam_names) > 0 & !all(model.tree.hp.names %in% names(kernel.hyper.params))) {
      stop(paste("Hyperparameter names not as expected. Expected: ",
                 paste(model.tree.hp.names, collapse="; "),
                 "Actual: ",
                 paste(names(kernel.hyper.params), collapse="; "),
                 sep=""))
    }
    kernel.hyper.params <- kernel.hyper.params[model.tree.hp.names]
    names(kernel.hyper.params) <- kernel$hyperparam_names

    if (return.grad) {
      # Set up a full array for all of the hyperparameters in the model
      num.training.points <- nrow(x)

      cov.mat.grad.array <- array(0, c(num.training.points,
                                       num.training.points,
                                       length(model.tree$all.hyper.params) + 1)
      )

      # Only update the hyperparameters which are used for this kernel.
      if (length(kernel.hyper.param.indices) > 0) {

        inst.grad.array <- cached_call(get_covariance_matrix_grad,
                                       kernel=kernel,
                                       x=x,
                                       sigma.n=0,
                                       hyper.params=kernel.hyper.params,
                                       cache=cache)[, , -1]


        cov.mat.grad.array[, , kernel.hyper.param.indices + 1] <- inst.grad.array
      }

      inst.cov.mats[[inst.id]] <- cov.mat.grad.array

    } else if (return.hess) {

      # Set up a full array for all of the hyperparameters in the model
      num.training.points <- nrow(x)

      cov.mat.hess.array <- array(0, c(num.training.points,
                                       num.training.points,
                                       length(model.tree$all.hyper.params) + 1,
                                       length(model.tree$all.hyper.params) + 1)
      )

      # Only update the hyperparameters which are used for this kernel.
      if (length(kernel.hyper.param.indices) > 0) {

        inst.hess.array <- cached_call(get_covariance_matrix_hess,
                                       kernel=kernel,
                                       x=x,
                                       sigma.n=0,
                                       hyper.params=kernel.hyper.params,
                                       cache=cache)[, , -1, -1]


        cov.mat.hess.array[, , kernel.hyper.param.indices + 1, kernel.hyper.param.indices + 1] <- inst.hess.array
      }

      inst.cov.mats[[inst.id]] <- cov.mat.hess.array


    } else {
      # Returning covariance matrix
      inst.cov.mats[[inst.id]] <- cached_call(get_covariance_matrix,
                                              kernel=kernel,
                                              x=x,
                                              sigma.n=0,
                                              hyper.params=kernel.hyper.params,
                                              cache=cache)
    }
  }

  return(inst.cov.mats)
}


#' Get the Covariance Matrix for a Kernel
#'
#' @param kernel A Kernel object
#' @param x A matrix of data to find the covariance matrix of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params Any additional parameters of \code{kernel} (internal use only)
#'
#' @return A covariance matrix
#' @export
get_covariance_matrix <- function(kernel, x, sigma.n, hyper.params, cache=NULL, ...) UseMethod("get_covariance_matrix")

#' Get the Covariance Matrix for a Kernel
#'
#' @param kernel A Kernel object
#' @param x A matrix of data to find the covariance matrix of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params Any additional parameters of \code{kernel} (internal use only)
#'
#' @return A covariance matrix
#' @export
get_covariance_matrix.Kernel <- function(kernel, x, sigma.n, hyper.params, cache=NULL) {
  k <- kernel$kernel
  additional.params <- kernel$additional_params
  return(get_covariance_matrix(kernel=k,
                               x=x,
                               sigma.n=sigma.n,
                               hyper.params=hyper.params,
                               additional.params=additional.params,
                               cache=cache))
}

#' Get the Covariance Matrix for a Kernel ModelTree
#'
#' @param kernel A ModelTree object
#' @param x A matrix of data to find the covariance matrix of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params an empty list (not used for ModelTree objects)
#'
#' @return A covariance matrix
#' @export

get_covariance_matrix.ModelTree <- function(kernel, x, sigma.n, hyper.params, additional.params=list(), cache=NULL) {
  model.tree <- kernel
  model.tree$all.hyper.params <- hyper.params

  inst.cov.mats <- get.inst.cov.mats(x, model.tree, return.type="cov", cache=cache)

  root.node <- find.root.node(model.tree)

  cov.mat <- eval.cov.mat.node(model.tree, inst.cov.mats, root.node)

  diag(cov.mat) <- diag(cov.mat) + sigma.n^2

  return(cov.mat)
}

#' Get the Covariance Matrix for a built-in kernel
#'
#' @param kernel A string
#' @param x A matrix of data to find the covariance matrix of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params Any additional parameters of \code{kernel} (internal use only)
#'
#' @return A covariance matrix
#' @export
get_covariance_matrix.character <- function(kernel, x, sigma.n, hyper.params, additional.params, cache=NULL) {
  if (length(hyper.params) == 0) {
    hyper.params <- numeric(0)
  }
  return(getCovarianceMatrixBuiltInCpp(x, kernel, sigma.n, hyper.params, additional.params))
}

#' Get the Covariance Matrix for a kernel function
#'
#' @param kernel A kernel function
#' @param x A matrix of data to find the covariance matrix of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params Any additional parameters of \code{kernel} (internal use only)
#'
#' @return A covariance matrix
#' @export
get_covariance_matrix.function <- function(kernel, x, sigma.n, hyper.params, additional.params, cache=NULL) {

  if (length(hyper.params) == 0) {
    hyper.params <- numeric(0)
  }

  num.training.points <- nrow(x)
  if (is.null(num.training.points)) {
    num.training.points <- length(x)
    x <- matrix(x, nrow=num.training.points)
  }

  cov.mat <- matrix(NA, nrow=num.training.points, ncol=num.training.points)
  for (sample.1 in 1:num.training.points) {
    for (sample.2 in 1:sample.1) {
      cov.mat[sample.1, sample.2] <- cov.mat[sample.2, sample.1] <-
        kernel(x[sample.1,], x[sample.2,], hyper.params, additional.params)
    }
  }

  diag(cov.mat) <- diag(cov.mat) + sigma.n^2

  return(cov.mat)
}

get_covariance_matrix_grad <- function(kernel, x, sigma.n, hyper.params, cache=NULL,  ...) UseMethod("get_covariance_matrix_grad")


#' Get the Covariance Matrix Grad for a Kernel
#'
#' @param kernel A Kernel object
#' @param x A matrix of data to find the covariance matrix grad of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params Any additional parameters of \code{kernel} (internal use only)
#'
#' @return A covariance matrix
#' @export
get_covariance_matrix_grad.Kernel <- function(kernel, x, sigma.n, hyper.params, cache=NULL) {
  k <- kernel$kernel
  k.grad <- kernel$grad
  additional.params <- kernel$additional_params
  return(get_covariance_matrix_grad(kernel=k,
                                    kernel.grad=k.grad,
                                    x=x,
                                    sigma.n=sigma.n,
                                    hyper.params=hyper.params,
                                    additional.params=additional.params,
                                    cache=cache))
}

#' Get the Covariance Matrix Grad for a Kernel ModelTree
#'
#' @param kernel A ModelTree object
#' @param kernel.grad Not used for ModelTree objects
#' @param x A matrix of data to find the covariance matrix grad of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params an empty list (not used for ModelTree objects)
#'
#' @return A covariance matrix
#' @export

get_covariance_matrix_grad.ModelTree <- function(kernel,
                                                 kernel.grad=NULL,
                                                 x,
                                                 sigma.n,
                                                 hyper.params,
                                                 additional.params=list(),
                                                 cache=NULL) {
  model.tree <- kernel
  model.tree$all.hyper.params <- hyper.params

  sigma.n.name <- "sigma.n"
  num.training.points <- nrow(x)

  inst.cov.mats <- get.inst.cov.mats(x, model.tree, return.type="cov", cache=cache)
  inst.cov.mat.grads <- get.inst.cov.mats(x, model.tree, return.type="grad", cache=cache)

  root.node <- find.root.node(model.tree)
  cov.mat.grad <- eval.cov.mat.grad.node(model.tree, inst.cov.mats, inst.cov.mat.grads, root.node)

  cov.mat.grad[, , 1] <- diag(2 * sigma.n, nrow=num.training.points)


  dimnames(cov.mat.grad) <- list(NULL,
                                 NULL,
                                 c(sigma.n.name,
                                   names(model.tree$all.hyper.params)
                                 )
  )

  return(cov.mat.grad)
}

#' Get the Covariance Matrix grad for a built-in kernel
#'
#' @param kernel A string
#' @param kernel.grad Not used for built in kernels
#' @param x A matrix of data to find the covariance matrix of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params Any additional parameters of \code{kernel} (internal use only)
#'
#' @return A covariance matrix
#' @export
get_covariance_matrix_grad.character <- function(kernel,
                                                 kernel.grad=NULL,
                                                 x,
                                                 sigma.n,
                                                 hyper.params,
                                                 additional.params,
                                                 cache=NULL) {
  # Make sure we pass a vector if there are no hyper params
  if (length(hyper.params) == 0) {
    hyper.params <- numeric(0)
  }

  # use the C++ function for built-in kernels
  cov.mat.grad.array <- getCovarianceMatrixGradArray(x,
                                                     kernel,
                                                     sigma.n,
                                                     hyper.params,
                                                     additional.params
  )

  dimnames(cov.mat.grad.array) <- list(NULL,
                                       NULL,
                                       c(sigma.n.name,
                                         names(hyper.params)
                                       )
  )

  return(cov.mat.grad.array)
}

#' Get the Covariance Matrix grad for a kernel function
#'
#' @param kernel A kernel function
#' @param kernel.grad A grad function for \code{kernel}
#' @param x A matrix of data to find the covariance matrix grad of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params Any additional parameters of \code{kernel} (internal use only)
#'
#' @return A covariance matrix
#' @export
get_covariance_matrix_grad.function <- function(kernel,
                                                kernel.grad,
                                                x,
                                                sigma.n,
                                                hyper.params,
                                                additional.params,
                                                cache=cache) {

  # Do the R thing
  num.training.points <- nrow(x)

  cov.mat.grad.array <- array(0, c(num.training.points,
                                   num.training.points,
                                   length(hyper.params) + 1)
  )


  cov.mat.grad.array[,,1] <- diag(2 * sigma.n, nrow=num.training.points)
  if (length(hyper.params) > 0) {
    for (sample.1 in 1:num.training.points) {
      for (sample.2 in 1:sample.1) {
        temp.grad <- kernel.grad(x[sample.1,], x[sample.2,], hyper.params, additional.params)
        cov.mat.grad.array[sample.1, sample.2, -1] <-
          cov.mat.grad.array[sample.2, sample.1, -1] <-
          temp.grad
      }
    }
  }

  dimnames(cov.mat.grad.array) <- list(NULL,
                                       NULL,
                                       c(sigma.n.name,
                                         names(hyper.params)
                                       )
  )

  return(cov.mat.grad.array)

}


# Hessian -----------------------------------------------------------------
get_covariance_matrix_hess <- function(kernel, x, sigma.n, hyper.params, cache=NULL,  ...) UseMethod("get_covariance_matrix_hess")


#' Get the Covariance Matrix Hessian for a Kernel
#'
#' @param kernel A Kernel object
#' @param x A matrix of data to find the covariance matrix Hessian of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params Any additional parameters of \code{kernel} (internal use only)
#'
#' @return A 4-d array containing the Hessian
#' @export
get_covariance_matrix_hess.Kernel <- function(kernel, x, sigma.n, hyper.params, cache=NULL) {
  k <- kernel$kernel
  k.hess <- kernel$hess
  additional.params <- kernel$additional_params
  return(get_covariance_matrix_hess(kernel=k,
                                    kernel.hess=k.hess,
                                    x=x,
                                    sigma.n=sigma.n,
                                    hyper.params=hyper.params,
                                    additional.params=additional.params,
                                    cache=cache))
}

#' Get the Covariance Matrix Hessian for a Kernel ModelTree
#'
#' @param kernel A ModelTree object
#' @param kernel.hess Not used for ModelTree objects
#' @param x A matrix of data to find the covariance matrix grad of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params an empty list (not used for ModelTree objects)
#'
#' @return A covariance matrix
#' @export

get_covariance_matrix_hess.ModelTree <- function(kernel,
                                                 kernel.hess=NULL,
                                                 x,
                                                 sigma.n,
                                                 hyper.params,
                                                 additional.params=list(),
                                                 cache=NULL) {
  model.tree <- kernel
  model.tree$all.hyper.params <- hyper.params

  sigma.n.name <- "sigma.n"
  num.training.points <- nrow(x)

  inst.cov.mats <- get.inst.cov.mats(x, model.tree, return.type="cov", cache=cache)
  inst.cov.mat.grads <- get.inst.cov.mats(x, model.tree, return.type="grad", cache=cache)
  inst.cov.mat.hess <- get.inst.cov.mats(x, model.tree, return.type="hess", cache=cache)

  root.node <- find.root.node(model.tree)
  cov.mat.hess <- eval.cov.mat.hess.node(model.tree,
                                         inst.cov.mats,
                                         inst.cov.mat.grads,
                                         inst.cov.mat.hess,
                                         root.node)

  cov.mat.hess[, , 1, 1] <- diag(2, nrow=num.training.points)


  dimnames(cov.mat.hess) <- list(NULL,
                                 NULL,
                                 c(sigma.n.name, names(model.tree$all.hyper.params)),
                                 c(sigma.n.name, names(model.tree$all.hyper.params))
  )

  return(cov.mat.hess)
}

#' Get the Covariance Matrix Hessian for a built-in kernel
#'
#' @param kernel A string
#' @param kernel.hess Not used for built in kernels
#' @param x A matrix of data to find the covariance matrix Hessian of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params Any additional parameters of \code{kernel} (internal use only)
#'
#' @return A covariance matrix
#' @export
get_covariance_matrix_hess.character <- function(kernel,
                                                 kernel.hess=NULL,
                                                 x,
                                                 sigma.n,
                                                 hyper.params,
                                                 additional.params,
                                                 cache=NULL) {
  # Make sure we pass a vector if there are no hyper params
  if (length(hyper.params) == 0) {
    hyper.params <- numeric(0)
  }

  # use the C++ function for built-in kernels
  cov.mat.hess.array <- getCovarianceMatrixHessianArray(x,
                                                  kernel,
                                                  sigma.n,
                                                  hyper.params,
                                                  additional.params
  )

  dimnames(cov.mat.hess.array) <- list(NULL,
                                       NULL,
                                       c(sigma.n.name, names(hyper.params)),
                                       c(sigma.n.name, names(hyper.params))
  )

  return(cov.mat.hess.array)
}

#' Get the Covariance Matrix Hessian for a kernel function
#'
#' @param kernel A kernel function
#' @param kernel.hess A Hessian function for \code{kernel}
#' @param x A matrix of data to find the covariance matrix grad of
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel}
#' @param additional.params Any additional parameters of \code{kernel} (internal use only)
#'
#' @return A covariance matrix
#' @export
get_covariance_matrix_hess.function <- function(kernel,
                                                kernel.hess,
                                                x,
                                                sigma.n,
                                                hyper.params,
                                                additional.params,
                                                cache=cache) {

  # Do the R thing
  num.training.points <- nrow(x)

  cov.mat.hess.array <- array(0, c(num.training.points,
                                   num.training.points,
                                   length(hyper.params) + 1,
                                   length(hyper.params) + 1)
  )


  cov.mat.grad.array[, , 1, 1] <- diag(2, nrow=num.training.points)
  if (length(hyper.params) > 0) {
    for (sample.1 in 1:num.training.points) {
      for (sample.2 in 1:sample.1) {
        temp.hess <- kernel.hess(x[sample.1,], x[sample.2,], hyper.params, additional.params)
        cov.mat.hess.array[sample.1, sample.2, -1, -1] <-
          cov.mat.hess.array[sample.2, sample.1, -1, -1] <-
          temp.grad
      }
    }
  }

  dimnames(cov.mat.hess.array) <- list(NULL,
                                       NULL,
                                       c(sigma.n.name, names(hyper.params)),
                                       c(sigma.n.name, names(hyper.params))
  )

  return(cov.mat.hess.array)

}

#' Get predictive K* matrix
#'
#' This function returns the K* matrix of covariances between the
#' training data and the data to predict
#'
#' @param kernel the kernel function
#' @param data.to.predict numeric matrix of data to predict
#' @param training.data the data used to train the Gaussian process
#' @param hyper.params the kernel's hyperparameters
#'
#' @return a numeric matrix
#' @export
get_kstar_matrix <- function(kernel,
                             data.to.predict,
                             training.data,
                             hyper.params,
                             ...) {
  UseMethod("get_kstar_matrix")
}

#' Get predictive K* matrix
#'
#' This function returns the K* matrix of covariances between the
#' training data and the data to predict for a Kernel object
#'
#' @param kernel the Kernel object
#' @param data.to.predict numeric matrix of data to predict
#' @param training.data the data used to train the Gaussian process
#' @param hyper.params the kernel's hyperparameters
#'
#' @return a numeric matrix
#' @export
get_kstar_matrix.Kernel <- function(kernel,
                             data.to.predict,
                             training.data,
                             hyper.params) {

  k <- kernel$kernel
  additional.params <- kernel$additional_params
  return(get_kstar_matrix(kernel=k,
                          data.to.predict=data.to.predict,
                          training.data=training.data,
                          hyper.params=hyper.params,
                          additional.params=additional.params))

}
#' Get predictive K* matrix
#'
#' This function returns the K* matrix of covariances between the
#' training data and the data to predict for an R function
#'
#' @param data.to.predict numeric matrix of data to predict
#' @param training.data the data used to train the Gaussian process
#' @param k the kernel function
#' @param hyper.params the kernel's hyperparameters
#' @param additional.params the kernel's additional parameters
#'
#' @return a numeric matrix
#' @export
get_kstar_matrix.function <- function(kernel,
                                      data.to.predict,
                                      training.data,
                                      hyper.params,
                                      additional.params) {

  num.data.points <- nrow(data.to.predict)
  num.training.points <- nrow(training.data)

  K.star <- matrix(NA, nrow=num.data.points, ncol=num.training.points)

  for (i in 1:num.data.points) {
    for (j in 1:num.training.points) {
      K.star[i, j] <- kernel(data.to.predict[i, ],
                             training.data[j, ],
                             hyper.params,
                             additional.params)
    }
  }

  return(K.star)
}

#' Get predictive K* matrix
#'
#' This function returns the K* matrix of covariances between the
#' training data and the data to predict for a built-in kernel
#'
#' @param kernel the kernel function
#' @param data.to.predict numeric matrix of data to predict
#' @param training.data the data used to train the Gaussian process
#' @param hyper.params the kernel's hyperparameters
#' @param additional.params the kernel's additional parameters
#'
#' @return a numeric matrix
#' @export
get_kstar_matrix.character <- function(kernel,
                                       data.to.predict,
                                       training.data,
                                       hyper.params,
                                       additional.params) {

  kernel_fun <- function(a,
                         b,
                         hyper.params,
                         additional.params) {
    callKernelByString(kernel, a, b, hyper.params, additional.params)
  }

  K.star <- get_kstar_matrix(kernel_fun,
                             data.to.predict,
                             training.data,
                             hyper.params,
                             additional.params)

  return(K.star)
}


#' K* matrix of a model tree node
#'
#' Recursively evaluates the K* matrix of a model tree from a given node.
#'
#' @param model.tree a model tree
#' @param inst.kstar.mats a named list of the K* matrices of the kernel instances in the model tree
#' @param node the node to evaluate
#'
#' @return a numeric matrix
eval.node.kstar.mat <- function(model.tree, inst.kstar.mats, node) {
  # More copypasta - maybe fix this too?
  node.entry <- model.tree$tree[node, ]
  if (!is.na(node.entry["operationID"])) {
    opID <- as.character(node.entry["operationID"])
    operation <- model.tree$operation.functions[[opID]]

    left.daughter <- as.numeric(node.entry["leftDaughter"])
    right.daughter <- as.numeric(node.entry["rightDaughter"])

    left.daughter.kstar.mat <- eval.node.kstar.mat(model.tree, inst.kstar.mats, left.daughter)
    right.daughter.kstar.mat <- eval.node.kstar.mat(model.tree, inst.kstar.mats, right.daughter)

    return(operation(left.daughter.kstar.mat, right.daughter.kstar.mat))
  } else if (!is.na(node.entry["kernelInstanceID"])) {
    inst.id <- as.character(node.entry["kernelInstanceID"])
    return(inst.kstar.mats[[inst.id]])
  } else {
    stop(paste("operationID and kernelInstanceID both NULL for node", node))
  }
}

#' K* Matrix of a Model Tree
#'
#' Returns the K* matrix between training data and data to predict for a given modelTree object.
#'
#' @param kernel a ModelTree object
#' @param data.to.predict numeric matrix of data to predict
#' @param training.data the data used to train the Gaussian process
#' @param hyper.params the ModelTree's hyperparameters
#' @param additional.params an empty list (not used for ModelTree objects)
#'
#' @return a numeric matrix
#' @export
get_kstar_matrix.ModelTree <- function(kernel,
                                       data.to.predict,
                                       training.data,
                                       hyper.params,
                                       additional.params=list()) {

  num.data.points <- nrow(data.to.predict)
  num.training.points <- nrow(training.data)

  #We re-use kernel later - beware!
  model.tree <- kernel

  model.tree$all.hyper.params <- hyper.params

  inst.kstar.mats <- list()

  # Copy pasta from the code for inst cov matrices.
  for (inst.id in rownames(model.tree$kernel.instances)) {
    kernel.class <- model.tree$kernel.instances[inst.id, "kernelClassName"]
    kernel <- model.tree$kernel.objects[[kernel.class]]
    kernel.hyper.param.indices <- model.tree$kernel.inst.hyper.param.indices[[inst.id]]
    kernel.hyper.params <- model.tree$all.hyper.params[kernel.hyper.param.indices]
    model.tree.hp.names <- paste(inst.id, kernel$hyperparam_names, sep="#")
    if (length(kernel$hyperparam_names) > 0 & !all(model.tree.hp.names %in% names(kernel.hyper.params))) {
      stop(paste("Hyperparameter names not as expected. Expected: ",
                 paste(model.tree.hp.names, collapse="; "),
                 "Actual: ",
                 paste(names(kernel.hyper.params), collapse="; "),
                 sep=""))
    }
    kernel.hyper.params <- kernel.hyper.params[model.tree.hp.names]
    names(kernel.hyper.params) <- kernel$hyperparam_names
    inst.kstar.mats[[inst.id]] <- get_kstar_matrix.Kernel(kernel,
                                                          data.to.predict,
                                                          training.data,
                                                          kernel.hyper.params)
  }

  root.node <- find.root.node(model.tree)

  kstar.mat <- eval.node.kstar.mat(model.tree, inst.kstar.mats, root.node)

  return(kstar.mat)
}
