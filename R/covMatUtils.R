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


#' Get the Covariance Matrix for a Kernel
#'
#' @param x A matrix of data to find the covariance matrix of
#' @param k A kernel function
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{k}
#' @param additional.params Any additional parameters of \code{k}
#'
#' @return A covariance matrix
get.covariance.matrix.kernel <- function(x, k, sigma.n, hyper.params, additional.params) {

  if (is.character(k)) {
    # Make sure we pass a vector if there are no hyper params
    if (length(hyper.params) == 0) {
      hyper.params <- numeric()
    }
    return(getCovarianceMatrixBuiltInCpp(x, k, sigma.n, hyper.params, additional.params))

  } else {

    num.training.points <- nrow(x)
    if (is.null(num.training.points)) {
      num.training.points <- length(x)
      x <- matrix(x, nrow=num.training.points)
    }

    cov.mat <- matrix(NA, nrow=num.training.points, ncol=num.training.points)
    for (sample.1 in 1:num.training.points) {
      for (sample.2 in 1:sample.1) {
        cov.mat[sample.1, sample.2] <- cov.mat[sample.2, sample.1] <-
          k(x[sample.1,], x[sample.2,], hyper.params, additional.params)
      }
    }
    diag(cov.mat) <- diag(cov.mat) + sigma.n^2
    return(cov.mat)

  }
}


#' Get the Covariance Matrix Grad for a Kernel
#'
#' @param x A matrix of data to find the covariance matrix of
#' @param kernel.grad A function returning the grad of a kernel (as a named vector)
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#' @param hyper.params The hyperparameters of \code{kernel.grad}
#' @param additional.params Any additional parameters of \code{kernel.grad}
#'
#' @return A covariance matrix grad (as an array)
get.covariance.matrix.grad.kernel <- function(x, kernel.grad, sigma.n, hyper.params, additional.params) {

  sigma.n.name <- "sigma.n"

  if (is.character(kernel.grad)) {
    # Make sure we pass a vector if there are no hyper params
    if (length(hyper.params) == 0) {
      hyper.params <- numeric()
    }
    # use the C++ function for built-in kernels
    cov.mat.grad.array <- getCovarianceMatrixGradArray(x,
                                                       kernel.grad,
                                                       sigma.n,
                                                       hyper.params,
                                                       additional.params
    )

  } else {
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
  }

  dimnames(cov.mat.grad.array) <- list(NULL,
                                       NULL,
                                       c(sigma.n.name,
                                         names(hyper.params)
                                       )
  )
  return(cov.mat.grad.array)
}


#' Get Kernel Instance Covariance Matrices or Grads
#'
#' Calculates the covariance matrices (or associated grads) for the kernel instances appearing in a model tree
#'
#' @param x A numeric matrix
#' @param model.tree A modelTree object
#' @param return.grad Whether to return the covariance matrices or the associated grads
#'
#' @return A named list containing covariance matrices (or grad arrays), with one entry for each kernel instance.
get.inst.cov.mats <- function(x, model.tree, return.grad=FALSE) {

  inst.cov.mats <- list()

  hyper.param.cache.check(model.tree)

  for (inst.id in rownames(model.tree$kernel.instances)) {
    if (return.grad) {
      # Set up a full array for all of the hyperparameters in the model
      num.training.points <- nrow(x)

      cov.mat.grad.array <- array(0, c(num.training.points,
                                       num.training.points,
                                       length(model.tree$all.hyper.params) + 1)
      )
    }
    kernel.class <- model.tree$kernel.instances[inst.id, "kernelClassName"]
    kernel <- model.tree$kernel.class.functions[[kernel.class]]
    kernel.hyper.param.indices <- model.tree$kernel.inst.hyper.param.indices[[inst.id]]
    kernel.hyper.param.names <- model.tree$kernel.class.hyper.param.names[[kernel.class]]
    kernel.hyper.params <- model.tree$all.hyper.params[kernel.hyper.param.indices]
    kernel.additional.params <- model.tree$kernel.class.additional.params[[kernel.class]]
    if (!is.null(kernel.hyper.param.names)) {
      names(kernel.hyper.params) <- kernel.hyper.param.names
    }
    if (return.grad) {

      # Only update the hyperparameters which are used for this kernel.
      if (length(kernel.hyper.param.indices) > 0) {

        # First check the cache
        if (inst.id %in% names(model.tree$cache[[icmg.cache.name]])) {
          inst.grad.array <- model.tree$cache[[icmg.cache.name]][[inst.id]]
        } else {
          inst.grad.array <-
            get.covariance.matrix.grad.kernel(x,
                                              kernel,
                                              0,
                                              kernel.hyper.params,
                                              kernel.additional.params)[, , -1]

          model.tree$cache[[icmg.cache.name]][[inst.id]] <- inst.grad.array
        }
        cov.mat.grad.array[,,kernel.hyper.param.indices + 1] <- inst.grad.array
      }

      inst.cov.mats[[inst.id]] <- cov.mat.grad.array

    } else {
      # Returning covariance matrix
      if (inst.id %in% names(model.tree$cache[[icm.cache.name]])) {
        inst.cov.mats[[inst.id]] <- model.tree$cache[[icm.cache.name]][[inst.id]]
      } else {
        inst.cov.mats[[inst.id]] <-
          get.covariance.matrix.kernel(x,
                                       kernel,
                                       0,
                                       kernel.hyper.params,
                                       kernel.additional.params)

        model.tree$cache[[icm.cache.name]][[inst.id]] <- inst.cov.mats[[inst.id]]
      }
    }
  }

  return(inst.cov.mats)
}


#' Covariance Matrix of a Model Tree
#'
#' Calculates the covariance matrix for a modelTree object
#'
#' @param x A numeric matrix
#' @param model.tree A modelTree object
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#'
#'
#' @return A covariance matrix
#' @export
#'
#' @examples
#' mt <- create.model.tree.builtin()
#' mt <- insert.kernel.instance(mt, 1, "SE", NULL, hyper.params=c(l=1))
#'
#' x <- rnorm(50)
#'
#' cov.mat <- get.covariance.matrix.model.tree(x, mt, sigma.n=0.1)
get.covariance.matrix.model.tree <- function(x, model.tree, sigma.n) {

  inst.cov.mats <- get.inst.cov.mats(x, model.tree, return.grad=FALSE)

  root.node <- find.root.node(model.tree)

  cov.mat <- eval.cov.mat.node(model.tree, inst.cov.mats, root.node)

  diag(cov.mat) <- diag(cov.mat) + sigma.n^2

  return(cov.mat)

}


#' Covariance Matrix Grad of a Model Tree
#'
#' Calculates the covariance matrix grad for a modelTree object
#'
#' @param x A numeric matrix
#' @param model.tree A modelTree object
#' @param sigma.n The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix
#'
#' @return An array containing the covariance matrix grad.
#' @export
#'
#' @examples
#' mt <- create.model.tree.builtin()
#' mt <- insert.kernel.instance(mt, 1, "SE", NULL, hyper.params=c(l=1))
#'
#' x <- rnorm(50)
#'
#' cov.mat <- get.covariance.matrix.model.tree(x, mt, sigma.n=0.1)
get.covariance.matrix.grad.model.tree <- function(x, model.tree, sigma.n) {
  sigma.n.name <- "sigma.n"
  num.training.points <- nrow(x)

  # Since we're doing this twice it would be good to build a cache for the
  # covariance matrices
  inst.cov.mats <- get.inst.cov.mats(x, model.tree, return.grad=FALSE)
  inst.cov.mat.grads <- get.inst.cov.mats(x, model.tree, return.grad=TRUE)

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

#' Get predictive K* matrix
#'
#' This function returns the K* matrix of covariances between the training data and the data to predict.
#'
#' @param data.to.predict numeric matrix of data to predict
#' @param training.data the data used to train the Gaussian process
#' @param k the kernel function
#' @param hyper.params the kernel's hyperparameters
#' @param additional.params the kernel's additional parameters
#'
#' @return a numeric matrix
#' @export
get.kstar.mat.kernel <- function(data.to.predict, training.data, k, hyper.params, additional.params) {

  num.data.points <- nrow(data.to.predict)
  num.training.points <- nrow(training.data)

  if (is.character(k)) {
    kernel <- function(a, b, hyper.params, additional.params) callKernelByString(k, a, b, hyper.params, additional.params)
  } else {
    kernel <- k
  }

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
#' @param data.to.predict a numeric matrix
#' @param training.data a numeric matrix
#' @param model.tree a modelTree object
#'
#' @return a numeric matrix
#' @export
get.kstar.mat.model.tree <- function(data.to.predict, training.data, model.tree) {

  num.data.points <- nrow(data.to.predict)
  num.training.points <- nrow(training.data)


  inst.kstar.mats <- list()

  # Copy pasta from the code for inst cov matrices.
  for (inst.id in rownames(model.tree$kernel.instances)) {
    kernel.class <- model.tree$kernel.instances[inst.id, "kernelClassName"]
    kernel <- model.tree$kernel.class.functions[[kernel.class]]
    kernel.hyper.param.indices <- model.tree$kernel.inst.hyper.param.indices[[inst.id]]
    kernel.hyper.params <- model.tree$all.hyper.params[kernel.hyper.param.indices]
    kernel.additional.params <- model.tree$kernel.class.additional.params[[kernel.class]]
    inst.kstar.mats[[inst.id]] <- get.kstar.mat.kernel(data.to.predict,
                                                       training.data,
                                                       kernel,
                                                       kernel.hyper.params,
                                                       kernel.additional.params)
  }

  root.node <- find.root.node(model.tree)

  kstar.mat <- eval.node.kstar.mat(model.tree, inst.kstar.mats, root.node)

  return(kstar.mat)
}
