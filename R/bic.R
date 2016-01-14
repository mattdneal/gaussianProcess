# Given a trained GP object, return the BIC value
bayesian.information.criterion <- function(gp.obj) {
  # Must be a trained GP
  log.L <- gp.obj$log.marginal.likelihood
  k <- length(gp.obj$model.tree$all.hyper.params) + 1
  n <- nrow(gp.obj$training.points)
  return(-2 * log.L + k * log(n))
}

# Take a model tree and some data, and build the model using BIC
bic.model.search <- function(x, y, base.model.tree, plot.gp=FALSE, reset.params=FALSE, max.new.nodes=3, ...) {
  model.bic.vec <- numeric()
  if (nrow(base.model.tree$tree) != 0) {
    # There is a base model defined, so first evaluate that.
    gp.obj <- create.gaussian.process(x, y, base.model.tree)
    result <- gp.obj$fit.hyperparams(sigma.n.init=NA, ...)

    model.bic.vec[as.character(gp.obj$model.tree)] <- bayesian.information.criterion(gp.obj)
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
