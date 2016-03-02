

#' Calculate BIC of a GP
#'
#' @param gp.obj A trained gaussianProcess object
#'
#' @return The BIC value for the input GP object.
#' @export
#' @seealso \code{\link{bic.model.search}} and \code{\link{gaussianProcess}}
#'
#' @importFrom cacheMan cached_call
#'
#' @examples
#' x <- rnorm(50)
#' y <- sin(1/(x^2 + 0.15))
#' mt <- create.model.tree.builtin()
#' mt <- insert.kernel.instance(mt, 1, "SE", NULL, hyper.params=c(l=1))
#' gp <- create.gaussian.process(x, y, mt)
#' gp$fit.hyperparams(NA)
bayesian.information.criterion <- function(gp.obj) {
  # Must be a trained GP
  log.L <- get.marginal.likelihood(gp=gp.obj,
                                   sigma.n=gp.obj$optimized.sigma.n,
                                   hyper.params=gp.obj$optimized.hyperparams)
  k <- length(gp.obj$optimized.hyperparams) + 1
  n <- nrow(gp.obj$training.points)
  return(-2 * log.L + k * log(n))
}

#' Calculate BIC of a GP (Effective DoF)
#'
#' Calculate the BIC using the effective degrees of freedom of the kernel,
#' rather than the number of free parameters
#'
#' @param gp.obj A trained gaussianProcess object
#'
#' @return The BIC value for the input GP object.
#' @export
#' @seealso \code{\link{bic.model.search}} and \code{\link{gaussianProcess}}
#'
#' @importFrom cacheMan cached_call
#'
#' @examples
#' x <- rnorm(50)
#' y <- sin(1/(x^2 + 0.15))
#' mt <- create.model.tree.builtin()
#' mt <- insert.kernel.instance(mt, 1, "SE", NULL, hyper.params=c(l=1))
#' gp <- create.gaussian.process(x, y, mt)
#' gp$fit.hyperparams(NA)
bayesian.information.criterion.EDoF <- function(gp.obj) {
  # Must be a trained GP
  log.L <- get.marginal.likelihood(gp=gp.obj,
                                   sigma.n=gp.obj$optimized.sigma.n,
                                   hyper.params=gp.obj$optimized.hyperparams)


  cov.mat.chol <- cached_call(get.cov.mat.chol,
                              kernel=gp.obj$kernel,
                              x=gp.obj$training.points,
                              sigma.n=gp.obj$optimized.sigma.n,
                              hyper.params=gp.obj$optimized.hyperparams,
                              cache=gp.obj$cache)

  cov.mat.inv <- cached_call(chol2inv,
                             cov.mat.chol,
                             cache=gp.obj$cache)

  cov.mat.noisefree <- cached_call(get_covariance_matrix,
                                   kernel=gp.obj$kernel,
                                   x=gp.obj$training.points,
                                   sigma.n=0,
                                   hyper.params=gp.obj$optimized.hyperparams,
                                   cache=gp.obj$cache)


  # Trace of matrix product is sum of pointwise product
  k <- sum(cov.mat.inv * cov.mat.noisefree)
  print(k)
  n <- nrow(gp.obj$training.points)
  return(-2 * log.L + k * log(n))
}


#' Calculate BIC of a GP (Effective DoF, eigenvalue decomposition)
#'
#' Calculate the BIC using the effective degrees of freedom of the kernel,
#' rather than the number of free parameters
#'
#' @param gp.obj A trained gaussianProcess object
#'
#' @return The BIC value for the input GP object.
#' @export
#' @seealso \code{\link{bic.model.search}} and \code{\link{gaussianProcess}}
#'
#' @importFrom cacheMan cached_call
#'
#' @examples
#' x <- rnorm(50)
#' y <- sin(1/(x^2 + 0.15))
#' mt <- create.model.tree.builtin()
#' mt <- insert.kernel.instance(mt, 1, "SE", NULL, hyper.params=c(l=1))
#' gp <- create.gaussian.process(x, y, mt)
#' gp$fit.hyperparams(NA)
bayesian.information.criterion.EDoF.eigen <- function(gp.obj) {
  # Must be a trained GP
  log.L <- get.marginal.likelihood(gp=gp.obj,
                                   sigma.n=gp.obj$optimized.sigma.n,
                                   hyper.params=gp.obj$optimized.hyperparams)
  
  cov.mat.noisefree <- cached_call(get_covariance_matrix,
                                   kernel=gp.obj$kernel,
                                   x=gp.obj$training.points,
                                   sigma.n=0,
                                   hyper.params=gp.obj$optimized.hyperparams,
                                   cache=gp.obj$cache)
  
  
  # Trace of matrix product is sum of pointwise product
  eigenvals <- eigen(cov.mat.noisefree, symmetric=TRUE, only.values=TRUE)$values
  k <- sum(eigenvals / (eigenvals + gp.obj$optimized.sigma.n^2))
  print(k)
  n <- nrow(gp.obj$training.points)
  return(-2 * log.L + k * log(n))
}


#' Perform Model Selection Using BIC
#'
#' @param x A matrix or data frame of predictors
#' @param y A numeric vector of responses
#' @param base.model.tree The model tree at which the search starts
#' @param plot.gp Whether to plot the gaussian processes encountered during the search. If ncol(x) > 1 this parameter is ignored and no plots are created.
#' @param reset.params By default the search starts with the optimum hyperparameter values found in the previous step in the search (with new hyperparameters fitted from random start points). Setting this to TRUE causes all hyperparameters to be randomly set at each step.
#' @param max.new.nodes The number of kernel instances to add to the base model.
#' @param revisit.kernels Whether to revisit previously-assessed kernels when deleting nodes from the current best candidate.
#' @param ... Additional parameters to be passed to gp.obj$fit.hyperparameters (see \code{\link{gaussianProcess}})
#'
#' @return The trained gaussianProcess object with the best BIC.
#' @export
#'
#' @importFrom cacheMan create_cache
#'
#' @examples
#' x <- rnorm(50)
#' y <- sin(1/(x^2 + 0.15))
#' mt <- create.model.tree.builtin()
#' gp <- bic.model.search(x, y, mt)
bic.model.search <- function(x,
                             y,
                             base.model.tree,
                             plot.gp=FALSE,
                             reset.params=FALSE,
                             max.new.nodes=3,
                             revisit.kernels=TRUE,
                             ...) {
  model.bic.vec <- numeric()

  evaluated.models <- character(0)

  cache <- create_cache()

  if (nrow(base.model.tree$tree) != 0) {
    # There is a base model defined, so first evaluate that.
    base.model.tree <- reduce.to.canonical.tree(base.model.tree)
    k <- create.kernel.object.from.model.tree(base.model.tree)
    gp.obj <- create.gaussian.process(x, y, k, cache)

    result <- fit.hyperparams(gp.obj, ...)
    gp.obj <- result$gp
    model.bic.vec[as.character(gp.obj$kernel$kernel)] <- bayesian.information.criterion(gp.obj)
    evaluated.models <- c(evaluated.models, as.character(gp.obj$kernel$kernel))
    print(paste("Base model:", as.character(gp.obj$kernel$kernel)))
    print(paste("Base model BIC:", model.bic.vec[as.character(gp.obj$kernel$kernel)]))

    best.bic <- model.bic.vec[as.character(gp.obj$kernel$kernel)]
    best.gp <- gp.obj
    best.model <- gp.obj$kernel$kernel
    if (plot.gp) plot(gp.obj)
  } else {
    best.bic <- Inf
    best.gp <- NULL
    best.model <- base.model.tree
  }

  bic.improved <- TRUE

  new.nodes <- 0

  while (bic.improved & new.nodes < max.new.nodes) {
    new.nodes <- new.nodes + 1

    if (reset.params & best.bic != Inf & any(!is.na(best.model$all.hyper.params))) {
      reset.vec <- c(FALSE, TRUE)
    } else {
      reset.vec <- FALSE
    }

    bic.improved <- FALSE
    next.models <- generate.next.models(best.model)

    if (!revisit.kernels) {
      duplicate.models <- numeric(0)

      for (i in seq_along(next.models)) {
        if (as.character(next.models[[i]]) %in% evaluated.models) {
          duplicate.models <- c(duplicate.models, i)
        }
      }
      if (length(duplicate.models) > 0) {
        next.models <- next.models[-duplicate.models]
      }
    }

    prev.best.hyper.params <- best.gp$optimized.hyperparams
    prev.best.sigma.n <- best.gp$optimized.sigma.n

    for (i in 1:length(next.models)) {
      for (reset in reset.vec) {
        current.model <- next.models[[i]]
        current.k <- create.kernel.object.from.model.tree(current.model)
        current.hyper.params <- rep(NA, length(current.k$hyperparam_names))
        names(current.hyper.params) <- current.k$hyperparam_names
        if (!reset) {
          current.hyper.params[names(prev.best.hyper.params)] <- prev.best.hyper.params
        }
        convcode <- 1
        counter <- 0
        while (convcode > 0 & counter < 3) {
          # This didn't converge - try again a few times, then give up and ignore this model
          counter <- counter + 1
          gp.obj <- create.gaussian.process(x, y, current.k, cache)
          result <- fit.hyperparams(gp.obj, hyper.params.init=current.hyper.params, ...)
          gp.obj <- result$gp
          convcode <- result$optimx.obj["convcode"]
        }

        if (result$optimx.obj["convcode"] == 0) {
          model.bic <- bayesian.information.criterion(gp.obj)
          model.bic.vec[paste(as.character(gp.obj$kernel$kernel), as.character(reset), sep=".")] <- model.bic
          evaluated.models <- c(evaluated.models, as.character(gp.obj$kernel$kernel))
          print(paste("Model:", as.character(gp.obj$kernel$kernel)))
          print(paste("Model BIC:", model.bic))
          print("Model hyperparams:")
          print(gp.obj$optimized.hyperparams)
          print(paste("sigma_n:", gp.obj$optimized.sigma.n))
          print(as.matrix(model.bic.vec[order(model.bic.vec)]))
          if (plot.gp) plot(gp.obj, main=as.character(gp.obj$kernel$kernel))
          if (model.bic < best.bic) {
            best.bic <- model.bic
            best.gp <- gp.obj
            best.model <- gp.obj$kernel$kernel
            bic.improved <- TRUE
          }
        } else {
          warnmess <- paste("Model", as.character(current.model), "failed to converge")
          warning(warnmess)
          print(warnmess)
        }
      }
    }
  }
  return(best.gp)
}
