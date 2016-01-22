
# By default generate up to 10^6, so that if we square the number we can
# still add on a small sigma_n and not have it round off
rplaw <- function(n, alpha, x_0=1, x_1=10^6) {
  alpha1 <- alpha + 1
  ((x_1^alpha1 - x_0^alpha1) * runif(n) + x_0^alpha1)^(1/alpha1)
}


#' Create a gaussianProcess object
#'
#' This function creates a gaussianProcess object along with associated functions for fitting hyperparameters.
#'
#' Intend to rework this to fit the R S3 framework - it currently encloses all its methods to allow caching in
#' optimx calls, but with a bit more experience I think this isn't actually necessary.
#'
#' @param x A matrix or data frame of predictor variables
#' @param y A numeric vector of response variables
#' @param model.tree The model tree defining the kernel function for the Gaussian process
#'
#' @return An untrained gaussianProcess object
#' @export
#'
#' @importFrom optimx optimx
#'
#' @examples
#' x <- rnorm(50)
#' y <- sin(1/(x^2 + 0.15))
#' mt <- create.model.tree.builtin()
#' mt <- insert.kernel.instance(mt, 1, "SE", NULL, hyper.params=c(l=1))
#' gp <- create.gaussian.process(x, y, mt)
#' gp$fit.hyperparams(NA)
#'
create.gaussian.process <- function(x, y, model.tree) {
  gp.obj <- new.env()
  class(gp.obj) <- "GaussianProcess"

  num.training.points <- nrow(x)
  if (is.null(num.training.points)) {
    num.training.points <- length(x)
    x <- matrix(x, nrow=num.training.points)
  }
  if (!is.numeric(num.training.points)) {
    stop("x does not appear to be a matrix or a vector")
  }

  gp.obj$training.points <- x
  gp.obj$num.training.points <- num.training.points
  gp.obj$training.point.values <- y
  gp.obj$model.tree <- clone.env(model.tree)

  gp.obj$saved.sigma.n <- NULL

  gp.obj$cov.mat <- NULL
  gp.obj$cov.mat.chol <- NULL
  gp.obj$cov.mat.L.inv <- NULL
  gp.obj$cov.mat.inv <- NULL
  gp.obj$alpha <- NULL
  gp.obj$cov.mat.det.log <- NULL

  gp.obj$log.marginal.likelihood <- NULL

  gp.obj$cov.mat.grad <- NULL
  gp.obj$log.marginal.likelihood.grad <- NULL

  # Train a GP with specific hyperparameters
  gp.obj$train <- function(sigma.n, hyper.params) {
    # Check if we've already done the calculations for these parameters, and if not then do them.
    recalculation.required <- FALSE

    if (is.null(model.tree$all.hyper.params)) {
      recalculation.required <- TRUE
    } else {
      # Check for NAs first
      if (any(is.na(model.tree$all.hyper.params)) | any(is.null(model.tree$all.hyper.params))) {
        recalculation.required <- TRUE
      } else if (any(hyper.params!=model.tree$all.hyper.params)) {
        recalculation.required <- TRUE
      }
    }

    if (is.null(saved.sigma.n)) {
      recalculation.required <- TRUE
    } else {
      if (saved.sigma.n != sigma.n) {
        recalculation.required <- TRUE
      }
    }

    if (recalculation.required) {


      #Save the hyper parameters
      saved.sigma.n <<- sigma.n
      names(hyper.params) <- names(model.tree$all.hyper.params)
      model.tree$all.hyper.params <<- hyper.params

      #Create the covariance matrix
      cov.mat <- get.covariance.matrix.model.tree(training.points, model.tree, sigma.n)

      #Store it in the enclosing environment
      cov.mat <<- cov.mat
      cov.mat.chol <<- chol(cov.mat)
      cov.mat.L.inv <<- solve(t(cov.mat.chol))
      cov.mat.inv <<- chol2inv(cov.mat.chol)
      alpha <<- cov.mat.inv %*% y
      cov.mat.det.log <<- 2 * sum(log(diag(cov.mat.chol)))
    }
  }
  environment(gp.obj$train) <- gp.obj

  # Find the marginal likelihood of the given hyperparameters
  gp.obj$get.marginal.likelihood <- function(sigma.n, hyper.params) {
    train(sigma.n, hyper.params)

    R <-  cov.mat.chol
    y <- training.point.values
    z <- forwardsolve(t(R), y)
    x <- backsolve(R, z)

    log.marginal.likelihood <<- as.numeric( - 1/2 * t(y) %*% x
                                            - 1/2 * cov.mat.det.log
                                            - num.training.points / 2
    )
    return(log.marginal.likelihood)
  }
  environment(gp.obj$get.marginal.likelihood) <- gp.obj

  # Find the grad of the marginal likelihood of the given hyperparameters
  gp.obj$get.marginal.likelihood.grad <- function(sigma.n, hyper.params) {
    train(sigma.n, hyper.params)

    sigma.n.name <- "sigma.n"

    cov.mat.grad.array <- get.covariance.matrix.grad.model.tree(training.points,
                                                                model.tree,
                                                                sigma.n
    )

    log.marginal.likelihood.grad <- numeric(length(hyper.params) + 1)
    names(log.marginal.likelihood.grad) <- c(sigma.n.name, names(hyper.params))

    # The * rather than %*% in what follows is not a mistake - tr(t(A)%*%B) is
    # the sum of the element-wise product of A and B.
    # Maybe replace the sum() with kahan's summation or something equivalent.
    # Can't find details of R's implementation of sum().
    temp.multiplicand <- t((alpha %*% t(alpha) - cov.mat.inv))

    temp.multiplicand <- array(rep(temp.multiplicand, length(hyper.params) + 1),
                               dim(cov.mat.grad.array)
    )

    #for (param.name in names(hyper.params)) {
    #  log.marginal.likelihood.grad[param.name] <- 1/2 * sum(temp.multiplicand * cov.mat.grad[[param.name]])
    #}

    #log.marginal.likelihood.grad[sigma.n.name] <- 1/2 * sum(temp.multiplicand * cov.mat.grad[[sigma.n.name]])

    log.marginal.likelihood.grad <- colSums(1/2 * temp.multiplicand * cov.mat.grad.array, dims=2)
    cov.mat.grad <<- cov.mat.grad
    log.marginal.likelihood.grad <<- log.marginal.likelihood.grad
    return(log.marginal.likelihood.grad)
  }
  environment(gp.obj$get.marginal.likelihood.grad) <- gp.obj

  #optimx wrapper for gp.obj$get.marginal.likelihood
  gp.obj$get.marginal.likelihood.optimx <- function(par) {
    -gp.obj$get.marginal.likelihood(sigma.n=par[1], hyper.params=par[-1])
  }
  environment(gp.obj$get.marginal.likelihood.optimx) <- gp.obj

  #optimx wrapper for gp.obj$get.marginal.likelihood.grad
  gp.obj$get.marginal.likelihood.grad.optimx <- function(par) {
    -gp.obj$get.marginal.likelihood.grad(sigma.n=par[1], hyper.params=par[-1])
  }
  environment(gp.obj$get.marginal.likelihood.grad.optimx) <- gp.obj

  gp.obj$fit.hyperparams <- function(sigma.n.init,
                                     random.init.gridsize=50,
                                     resample.top.inits=TRUE,
                                     resample.from=NULL,
                                     num.resamples=NULL,
                                     additional.optimx.runs=5,
                                     additional.par.perc.diff=0.1,
                                     random.init.scale=1,
                                     abs.min.sigma.n=0.001,
                                     abs.min.hyper.params=0,
                                     optimx.method="L-BFGS-B", max.iterations=10000,
                                     optimx.starttests=FALSE, optimx.trace=0,
                                     verbose=FALSE) {
    hyper.params.init <- model.tree$all.hyper.params

    # generate random.init.gridsize random pars and find the best
    best.par <- NULL
    best.log.likelihood <- -Inf

    par.mat <- matrix(0, nrow=random.init.gridsize, ncol=length(hyper.params.init) + 1)
    colnames(par.mat) <- c("sigma.n", names(hyper.params.init))
    log.likelihood.vec <- rep(-Inf, random.init.gridsize)

    for (i in 1:random.init.gridsize) {
      if (is.na(sigma.n.init)) {
        sigma.n.temp <- abs(rcauchy(n=1, scale=random.init.scale)) + abs.min.sigma.n
      } else {
        sigma.n.temp <- sigma.n.init
      }
      if (length(hyper.params.init) > 0) {
        par <- c(sigma.n.temp, sapply(hyper.params.init, function(x) {
                                                          if (is.na(x)) {
                                                            randnum <- rcauchy(n=1, scale=random.init.scale)
                                                            return(sign(randnum) * (abs(randnum) + abs.min.hyper.params))
                                                          } else {
                                                            return(x)
                                                          }
                                                        }
                                     )
                )
      } else {
        par <- sigma.n.temp
      }
      names(par) <- c("sigma.n", names(hyper.params.init))

      # We're going to pick some crazy hyperparameters because we're randomly generating
      # them, and some of them are going to break the kernels (so far mainly by causing round-off
      # errors between the kernel and sigma_n). To avoid this, we're going to convert errors into
      # warnings, pretend the bad parameters never happened, and carry on with our day.
      tryCatch({
        log.likelihood.temp <- get.marginal.likelihood(par[1], par[-1])

        par.mat[i, ] <- par
        log.likelihood.vec[i] <- log.likelihood.temp

        if (log.likelihood.temp > best.log.likelihood | is.null(best.log.likelihood)) {
          best.log.likelihood <- log.likelihood.temp
          best.par <- par
        }
      },
      error=function(e) {
        warning(paste("Bad random hyperparams encountered. Error:", e))
      }
      )
    }
    if (verbose) print("Best pars found:")
    if (verbose) print(best.par)
    if (verbose) print(paste("Log likelihood:", best.log.likelihood))

    if (resample.top.inits) {
      if (is.null(resample.from)) {
        resample.from <- ceiling(random.init.gridsize * 0.1)
      }
      if (is.null(num.resamples)) {
        num.resamples <- random.init.gridsize
      }
      if (verbose) print(paste("Resampling based on top", resample.from, "candidate starting hyperparameters."))
      top.n.pars <- par.mat[order(log.likelihood.vec, decreasing=T)[1:resample.from], , drop=FALSE]

      if (verbose) print("Summary of top n hyperparam starting values:")
      if (verbose) print(summary(top.n.pars))

      resampled.par.mat <- matrix(0, nrow=num.resamples, ncol=length(hyper.params.init) + 1)
      resampled.log.likelihood.vec <- rep(-Inf, num.resamples)

      for (col in 1:ncol(top.n.pars)) {
        # Find a suitable bandwidth and then use this to sample from the
        # distribution of the top n for this par.
        tryCatch({
          bw <- bw.SJ(top.n.pars[, col])
        },
        error=function(e) {
          warning(e)
          bw <- bw.nrd0(top.n.pars[, col])
        })
        resampled.par.mat[, col] <- rnorm(num.resamples,
                                          sample(top.n.pars[, col],
                                                 size = num.resamples,
                                                 replace = TRUE),
                                          bw)
      }

      if (verbose) print(paste("Taking", num.resamples, "resamples of top n starting values."))
      for (i in 1:num.resamples) {
        par <- resampled.par.mat[i, ]
        names(par) <- c("sigma.n", names(hyper.params.init))
        # We're going to pick some crazy hyperparameters because we're randomly generating
        # them, and some of them are going to break the kernels (so far mainly by causing round-off
        # errors between the kernel and sigma_n). To avoid this, we're going to convert errors into
        # warnings, pretend the bad parameters never happened, and carry on with our day.
        tryCatch({
          log.likelihood.temp <- get.marginal.likelihood(par[1], par[-1])

          resampled.log.likelihood.vec[i] <- log.likelihood.temp

          if (log.likelihood.temp > best.log.likelihood | is.null(best.log.likelihood)) {
            best.log.likelihood <- log.likelihood.temp
            best.par <- par
          }
        },
        error=function(e) {
          warning(paste("Bad random hyperparams encountered. Error:", e))
        }
        )
      }
      if (verbose) print("Best pars found:")
      if (verbose) print(best.par)
      if (verbose) print(paste("Log likelihood:", best.log.likelihood))
      par.mat <- rbind(par.mat, resampled.par.mat)
      log.likelihood.vec <- c(log.likelihood.vec, resampled.log.likelihood.vec)
    }
    lower <- -Inf
    if (!is.null(abs.min.sigma.n) | abs.min.sigma.n != 0) {
      lower <- rep(-Inf, length(best.par))
      lower[1] <- abs.min.sigma.n
    }
    optimx.obj <- optimx(best.par, get.marginal.likelihood.optimx,
                         gr=get.marginal.likelihood.grad.optimx,
                         method=optimx.method,
                         itnmax=max.iterations,
                         lower=lower,
                         control=list(starttests=optimx.starttests,
                                      trace=optimx.trace)
    )

    if (optimx.obj["convcode"] > 0) {
      warning("Did not converge!")
      best.log.likelihood <- -Inf
    } else {
      opt.pars <- as.numeric(optimx.obj)[1:(length(hyper.params.init) + 1)]
      names(opt.pars) <- c("sigma.n", names(hyper.params.init))
      train(opt.pars[1], opt.pars[-1])
      best.log.likelihood <- gp.obj$log.marginal.likelihood
    }

    tested.par.indices <- c(1)
    candidate.par.index <- 1
    par.order <- order(log.likelihood.vec, decreasing=T)

    for (i in seq_len(additional.optimx.runs)) {
      # Take the next set of pars which is sufficiently different from the others tried
      par.diff <- 0
      while (candidate.par.index < nrow(par.mat) & par.diff < additional.par.perc.diff) {
        candidate.par.index <- candidate.par.index + 1
        new.par <- par.mat[par.order[candidate.par.index], ]
        par.diff <- min(apply(par.mat[tested.par.indices, , drop=FALSE],
                              1,
                              function(row) sqrt(sum((new.par - row)^2)) / sqrt(sum((row)^2))
        ))
      }

      if (candidate.par.index >= nrow(par.mat) | log.likelihood.vec[candidate.par.index] == -Inf) {
        # We've not got a sufficiently different set of pars
        break
      }

      tested.par.indices <- c(tested.par.indices, candidate.par.index)

      new.optimx.obj <- optimx(new.par, get.marginal.likelihood.optimx,
                               gr=get.marginal.likelihood.grad.optimx,
                               method=optimx.method,
                               itnmax=max.iterations,
                               lower=lower,
                               control=list(starttests=optimx.starttests,
                                            trace=optimx.trace)
      )

      if (new.optimx.obj["convcode"] > 0) {
        warning("Did not converge!")
      } else {
        trained.pars <- as.numeric(new.optimx.obj)[1:(length(hyper.params.init) + 1)]
        names(trained.pars) <- c("sigma.n", names(hyper.params.init))
        train(trained.pars[1], trained.pars[-1])
        if (gp.obj$log.marginal.likelihood > best.log.likelihood) {
          best.log.likelihood <- gp.obj$log.marginal.likelihood
          opt.pars <- trained.pars
          optimx.obj <- new.optimx.obj
        }
      }

    }
    if (best.log.likelihood != -Inf) {
      train(opt.pars[1], opt.pars[-1])
      if (verbose) print(paste("Final.pars after optimisation:"))
      if (verbose) print(opt.pars)
      if (verbose) print(paste("Log likelihood for best starting par optimisation:", best.log.likelihood))
    }
    return(optimx.obj)
  }
  environment(gp.obj$fit.hyperparams) <- gp.obj

  return(gp.obj)
}


#' Create Predictions Using a Gaussian Process
#'
#' Returns the posterior mean and variance of a set of data points under a given Gaussian process
#'
#' @param gp.obj A trained gaussianProcess object
#' @param data Data to predict
#'
#' @return a list containing named elements:
#' \itemize{
#'   \item \code{mean} - the predicted mean value for each data point in \code{data}.
#'   \item \code{var} - the variance about the predicted mean for each data point in \code{data}.
#' }
#' @export
#'
#' @examples
#' x <- rnorm(50)
#' y <- sin(1/(x^2 + 0.15))
#' mt <- create.model.tree.builtin()
#' mt <- insert.kernel.instance(mt, 1, "SE", NULL, hyper.params=c(l=1))
#' gp <- create.gaussian.process(x, y, mt)
#' gp$fit.hyperparams(NA)
#'
#' x1 <- rnorm(50)
#' y1.predicted <- predict(gp, x1)$mean
#'
predict.GaussianProcess <- function(gp.obj, data) {
  if (!is.matrix(data)) {
    data <- matrix(data, nrow=length(data))
  }

  num.data.points <- nrow(data)
  kernel.func.list <- list()
  K.star.inst.list <- list()

  K.star <- get.kstar.mat.model.tree(data, gp.obj$training.points, gp.obj$model.tree)

  pred.out <- list()
  pred.out[["mean"]] <- K.star %*% gp.obj$alpha

  pred.out[["var"]] <- matrix(NA, nrow=num.data.points, ncol=1)

  for (i in seq(num.data.points)) {
    v <- gp.obj$cov.mat.L.inv %*% t(K.star[i, , drop=F])
    pred.out$var[i, 1] <- get.covariance.matrix.model.tree(data[i, , drop=FALSE], gp.obj$model.tree, 0) - t(v) %*% v
  }

  return(pred.out)
}

#' Plot a Gaussian Process
#'
#' Plots the posterior mean and first three standard deviations of a one-dimensional Gaussian process.
#'
#' @param gp.obj A gaussianProcess object
#' @param xlim The limits of the x-axis in the plot. Defaults to the limits of the training data plus one-quarter of the range at either end.
#' @param num.points The number of data points to calculate for the plot
#' @param ... Additional parameters to pass to \code{plot}
#'
#' @return NULL
#' @export
plot.GaussianProcess <- function(gp.obj, xlim=NULL, num.points=1000, ...) {
  if (is.null(xlim)) {
    x.min <- min(gp.obj$training.points)
    x.max <- max(gp.obj$training.points)

    x.range <- x.max - x.min

    x.min <- x.min - x.range / 4
    x.max <- x.max + x.range / 4
  } else {
    x.min = xlim[1]
    x.max = xlim[2]
  }

  x.values <- seq(from=x.min, to=x.max, length.out=num.points)
  preds <- predict(gp.obj, x.values)

  lower.bounds <- preds$mean - 3*sqrt(preds$var)
  upper.bounds <- preds$mean + 3*sqrt(preds$var)

  lower.2sd <- preds$mean - 2*sqrt(preds$var)
  upper.2sd <- preds$mean + 2*sqrt(preds$var)

  lower.1sd <- preds$mean - 1*sqrt(preds$var)
  upper.1sd <- preds$mean + 1*sqrt(preds$var)


  lower.bound <- min(lower.bounds, gp.obj$training.point.values)
  upper.bound <- max(upper.bounds, gp.obj$training.point.values)

  if (is.na(lower.bound) | is.na(upper.bound)) {
    warning("Lower or Upper bound is NA")
  } else {
    plot(x.values, preds$mean, type="n", ylim=c(lower.bound, upper.bound), ...)


    polygon(c(rev(x.values), x.values), c(rev(lower.bounds), upper.bounds), col = 'grey80', border = NA)
    polygon(c(rev(x.values), x.values), c(rev(lower.2sd), upper.2sd), col = 'grey70', border = NA)
    polygon(c(rev(x.values), x.values), c(rev(lower.1sd), upper.1sd), col = 'grey60', border = NA)
    lines(x.values, preds$mean)

    points(gp.obj$training.points, gp.obj$training.point.values)
  }

}

#' Sample a function from a Gaussian Process Prior
#'
#' @param x Points to sample at
#' @param kernel Kernel of the Gaussian Process
#' @param sigma.n Standard deviation of the noise
#' @param hyper.params Hyperparameters of the kernel
#' @param num.functions The number of different functions to generate
#'
#' @return A numeric vector of values at the sampled points
#' @export
#'
#' @importFrom mnormt rmnorm
sample.functions.from.kernel <- function(x, kernel, sigma.n, hyper.params, num.functions=1) {
  cov.mat <- get.covariance.matrix.kernel(x, kernel, sigma.n, hyper.params)
  y <- rmnorm(n=num.functions, mean=rep(0, nrow(cov.mat)), cov.mat)
  return(y)
}


#' Sample a function from a Gaussian Process Prior
#'
#' @param x Points to sample at
#' @param model.tree Model tree which defines the kernel of the Gaussian Process
#' @param sigma.n Standard deviation of the noise
#' @param num.functions The number of different functions to generate
#'
#' @return A numeric vector of values at the sampled points
#' @export
#'
#' @importFrom mnormt rmnorm
sample.functions.from.model.tree <- function(x, model.tree, sigma.n, num.functions=1) {
  if (!is.matrix(x)) {
    x <- matrix(x)
  }
  cov.mat <- get.covariance.matrix.model.tree(x, model.tree, sigma.n)
  y <- rmnorm(n=num.functions, mean=rep(0, nrow(cov.mat)), cov.mat)
  return(y)
}



