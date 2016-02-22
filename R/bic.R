

#' Calculate BIC of a GP
#'
#' @param gp.obj A trained gaussianProcess object
#'
#' @return The BIC value for the input GP object.
#' @export
#' @seealso \code{\link{bic.model.search}} and \code{\link{gaussianProcess}}
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
  log.L <- gp.obj$log.marginal.likelihood
  k <- length(gp.obj$model.tree$all.hyper.params) + 1
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

  if (nrow(base.model.tree$tree) != 0) {
    # There is a base model defined, so first evaluate that.
    order.tree(base.model.tree)
    gp.obj <- create.gaussian.process(x, y, base.model.tree)
    result <- gp.obj$fit.hyperparams(sigma.n.init=NA, ...)

    model.bic.vec[as.character(gp.obj$model.tree)] <- bayesian.information.criterion(gp.obj)
    evaluated.models <- c(evaluated.models, as.character(gp.obj$model.tree))
    print(paste("Base model:", as.character(gp.obj$model.tree)))
    print(paste("Base model BIC:", model.bic.vec[as.character(gp.obj$model.tree)]))

    best.bic <- model.bic.vec[as.character(gp.obj$model.tree)]
    best.gp <- gp.obj
    best.model <- gp.obj$model.tree
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

    duplicate.models <- numeric(0)

    for (i in seq_along(next.models)) {
      if (as.character(next.models[[i]]) %in% evaluated.models) {
        duplicate.models <- c(duplicate.models, i)
      }
    }
    if (length(duplicate.models) > 0) {
      next.models <- next.models[-duplicate.models]
    }

    for (i in 1:length(next.models)) {
      for (reset in reset.vec) {
        current.model <- clone.env(next.models[[i]])
        if (reset) {
          current.model$all.hyper.params[] <- NA
        }
        gp.obj <- create.gaussian.process(x, y, current.model)
        result <- gp.obj$fit.hyperparams(sigma.n.init=NA, ...)

        counter <- 0
        while (result["convcode"] > 0 & counter < 3) {
          # This didn't converge - try again a few times, then give up and ignore this model
          counter <- counter + 1
          gp.obj <- create.gaussian.process(x, y, current.model)
          result <- gp.obj$fit.hyperparams(sigma.n.init=NA, ...)
        }

        if (result["convcode"] == 0) {
          model.bic <- bayesian.information.criterion(gp.obj)
          model.bic.vec[paste(as.character(gp.obj$model.tree), as.character(reset), sep=".")] <- model.bic
          evaluated.models <- c(evaluated.models, as.character(gp.obj$model.tree))
          print(paste("Model:", as.character(gp.obj$model.tree)))
          print(paste("Model BIC:", model.bic))
          print("Model hyperparams:")
          print(gp.obj$model.tree$all.hyper.params)
          print(paste("sigma_n:", gp.obj$saved.sigma.n))
          print(as.matrix(model.bic.vec[order(model.bic.vec)]))
          if (plot.gp) plot(gp.obj, main=as.character(gp.obj$model.tree))
          if (model.bic < best.bic) {
            best.bic <- model.bic
            best.gp <- gp.obj
            best.model <- gp.obj$model.tree
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
