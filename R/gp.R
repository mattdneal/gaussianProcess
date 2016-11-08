GaussianProcess_class_name <- "GaussianProcess"
sigma.n.name <- "sigma.n"

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
#' @param kernel a Kernel object specifying the kernel for the Gaussian process
#' @param cache a Cache object. If NULL, a cache is created
#'
#' @return An untrained gaussianProcess object
#' @export
#'
#' @importFrom optimx optimx
#' @importFrom cacheMan create_cache hash_object
#'
#' @examples
#' x <- rnorm(50)
#' y <- sin(1/(x^2 + 0.15))
#' mt <- create.model.tree.builtin()
#' mt <- insert.kernel.instance(mt, 1, "squaredExponential", NULL, hyper.params=c(l=NULL))
#' k <- create.kernel.object.from.model.tree(mt)
#' gp <- create.gaussian.process(x, y, k)
#' gp <- fit.hyperparams(gp)
#'
create.gaussian.process <- function(x, y, kernel, cache=NULL) {
  gp.obj <- list()
  class(gp.obj) <- GaussianProcess_class_name

  num.training.points <- nrow(x)
  if (is.null(num.training.points)) {
    num.training.points <- length(x)
    x <- matrix(x, nrow=num.training.points)
  }
  if (!is.numeric(num.training.points)) {
    stop("x does not appear to be a matrix or a vector")
  }

  if (class(kernel) != Kernel_class_name) {
    stop("invalid kernel specified")
  }

  gp.obj$training.points <- x
  gp.obj$num.training.points <- num.training.points
  gp.obj$training.point.values <- y
  gp.obj$kernel <- kernel

  if (is.null(cache)) {
    gp.obj$cache <- create_cache()
  } else {
    gp.obj$cache <- cache
  }

  attr(gp.obj$training.points, cacheMan:::hash_attr) <-
    hash_object(gp.obj$training.points, gp.obj$cache)

  attr(gp.obj$training.point.values, cacheMan:::hash_attr) <-
    hash_object(gp.obj$training.point.values, gp.obj$cache)

  gp.obj$kernel <- set.kernel.hash(gp.obj$kernel, cache)

  return(gp.obj)
}

check.hyperparams <- function(gp, hyper.params) {
  if (class(gp) != GaussianProcess_class_name) {
    stop("Invalid Gaussian process supplied.")
  }
  k.hp.names <- gp$kernel$hyperparam_names
  if (length(hyper.params) != length(k.hp.names)) {
    stop("Number of hyperparameters supplied does not match kernel of GP")
  }
  if (!is.null(names(hyper.params))) {
    if (any(!names(hyper.params) %in% k.hp.names)) {
      stop("Named vector of hyperparameters supplied does not match hyperparameter names of kernel")
    }
    hyper.params <- hyper.params[k.hp.names]
  } else {
    warning("Unnamed vector of hyperparams supplied.")
    names(hyper.params) <- k.hp.names
  }
  return(hyper.params)
}

#' Get the Cholesky decomposition of a covariance matrix
#'
#' @return the Cholesky decomposition of a covariance matrix
#'
#' @importFrom cacheMan cached_call
get.cov.mat.chol <- function(kernel,
                             x,
                             sigma.n,
                             hyper.params,
                             cache
) {

  cov.mat <- cached_call(get_covariance_matrix,
                         kernel=kernel,
                         x=x,
                         sigma.n=sigma.n,
                         hyper.params=hyper.params,
                         cache=cache)

  cov.mat.chol <- cached_call(chol,
                              x=cov.mat,
                              cache=cache)

  return(cov.mat.chol)
}

#' Get alpha for a GP
#'
#' Internal shenanigans. This documentation only exists so I can tell Roxygen
#' what to import for it.
#'
#' @param gp a GaussianProcess object
#'
#' @return a vector
#'
#' @importFrom cacheMan cached_call
get.alpha <- function(kernel,
                      training.points,
                      training.point.values,
                      sigma.n,
                      hyper.params,
                      cache) {

  cov.mat.chol <- cached_call(get.cov.mat.chol,
                              kernel=kernel,
                              x=training.points,
                              sigma.n=sigma.n,
                              hyper.params=hyper.params,
                              cache=cache
                              )

  cov.mat.inv <- cached_call(chol2inv,
                             cov.mat.chol,
                             cache=cache)

  y <- training.point.values

  alpha <- cov.mat.inv %*% y

  return(alpha)
}

#' Get Log Marginal Likelihood of a GP
#'
#' @param gp a GaussianProcess object
#' @param sigma.n a numeric value for sigma.n
#' @param hyper.params a numeric vector of hyperparameters
#'
#' @return the log marginal likelihood of the GP at the given hyperparameters
#'
#' @importFrom cacheMan cached_call
get.marginal.likelihood <- function(gp, sigma.n, hyper.params) {
  hyper.params <- check.hyperparams(gp, hyper.params)

  cov.mat.chol <- cached_call(get.cov.mat.chol,
                              kernel=gp$kernel,
                              x=gp$training.points,
                              sigma.n=sigma.n,
                              hyper.params=hyper.params,
                              cache=gp$cache)

  cov.mat.det.log <- 2 * sum(log(diag(cov.mat.chol)))

  R <- cov.mat.chol

  y <- gp$training.point.values

  z <- cached_call(forwardsolve,
                   t(R), y,
                   cache=gp$cache)

  x <- cached_call(backsolve,
                   R, z,
                   cache=gp$cache)

  log.marginal.likelihood <- as.numeric( - 1/2 * t(y) %*% x
                                         - 1/2 * cov.mat.det.log
                                         - nrow(gp$training.points)
                                           * log(2 * pi) / 2)

  return(log.marginal.likelihood)
}

#' Get the Grad of the Log Marginal Likelihood of a GP
#'
#' Returns the grad of the log marginal likelihood of a GP w.r.t. the kernel's
#' hyperparameters.
#'
#' @param gp a GaussianProcess object
#' @param sigma.n a numeric value for sigma.n
#' @param hyper.params a numeric vector of hyperparameters
#'
#' @return
#'
#' @importFrom cacheMan cached_call
#'
#' @examples
get.marginal.likelihood.grad <- function(gp, sigma.n, hyper.params) {
  hyper.params <- check.hyperparams(gp, hyper.params)

  cov.mat.grad.array <- cached_call(get_covariance_matrix_grad,
                                    kernel=gp$kernel,
                                    x=gp$training.points,
                                    sigma.n=sigma.n,
                                    hyper.params=hyper.params,
                                    cache=gp$cache)

  log.marginal.likelihood.grad <- numeric(length(hyper.params) + 1)
  names(log.marginal.likelihood.grad) <- c(sigma.n.name, names(hyper.params))

  cov.mat.chol <- cached_call(get.cov.mat.chol,
                              kernel=gp$kernel,
                              x=gp$training.points,
                              sigma.n=sigma.n,
                              hyper.params=hyper.params,
                              cache=gp$cache)

  cov.mat.inv <- cached_call(chol2inv,
                             cov.mat.chol,
                             cache=gp$cache)

  alpha <- cached_call(get.alpha,
                       kernel=gp$kernel,
                       training.points=gp$training.points,
                       training.point.values=gp$training.point.values,
                       sigma.n=sigma.n,
                       hyper.params=hyper.params,
                       cache=gp$cache)

  # The * rather than %*% in what follows is not a mistake - tr(t(A)%*%B) is
  # the sum of the element-wise product of A and B.
  # Maybe replace the sum() with kahan's summation or something equivalent.
  # Can't find details of R's implementation of sum().
  temp.multiplicand <- t((alpha %*% t(alpha) - cov.mat.inv))

  temp.multiplicand <- array(rep(temp.multiplicand, length(hyper.params) + 1),
                             dim(cov.mat.grad.array))

  log.marginal.likelihood.grad <- colSums(1/2 * temp.multiplicand * cov.mat.grad.array, dims=2)

  return(log.marginal.likelihood.grad)
}


#' Get the Hessian of the Log Marginal Likelihood of a GP
#'
#' Returns the Hessian matrix of the log marginal likelihood of a GP w.r.t. the kernel's
#' hyperparameters.
#'
#' @param gp a GaussianProcess object
#' @param sigma.n a numeric value for sigma.n
#' @param hyper.params a numeric vector of hyperparameters
#'
#' @return
#'
#' @importFrom cacheMan cached_call
#'
#' @examples
get.marginal.likelihood.hessian <- function(gp, sigma.n, hyper.params) {
  hyper.params <- check.hyperparams(gp, hyper.params)

  num.hps <- length(hyper.params) + 1

  cov.mat.grad.array <- cached_call(get_covariance_matrix_grad,
                                    kernel=gp$kernel,
                                    x=gp$training.points,
                                    sigma.n=sigma.n,
                                    hyper.params=hyper.params,
                                    cache=gp$cache)

  cov.mat.hess.array <- cached_call(get_covariance_matrix_hess,
                                    kernel=gp$kernel,
                                    x=gp$training.points,
                                    sigma.n=sigma.n,
                                    hyper.params=hyper.params,
                                    cache=gp$cache)

  log.marginal.likelihood.hess <- matrix(0,
                                         nrow=num.hps,
                                         ncol=num.hps)

  colnames(log.marginal.likelihood.hess) <- c(sigma.n.name, names(hyper.params))
  rownames(log.marginal.likelihood.hess) <- c(sigma.n.name, names(hyper.params))

  cov.mat.chol <- cached_call(get.cov.mat.chol,
                              kernel=gp$kernel,
                              x=gp$training.points,
                              sigma.n=sigma.n,
                              hyper.params=hyper.params,
                              cache=gp$cache)

  cov.mat.inv <- cached_call(chol2inv,
                             cov.mat.chol,
                             cache=gp$cache)

  alpha <- cached_call(get.alpha,
                       kernel=gp$kernel,
                       training.points=gp$training.points,
                       training.point.values=gp$training.point.values,
                       sigma.n=sigma.n,
                       hyper.params=hyper.params,
                       cache=gp$cache)


  grad_k_alpha <- matrix(0, nrow=length(alpha), ncol=num.hps)
  grad.k.k.inv <- cov.mat.grad.array
  for (i in 1:num.hps) {
    grad_k_alpha[, i] <- cov.mat.grad.array[, , i] %*% alpha
    grad.k.k.inv[, , i] <- cov.mat.grad.array[, , i] %*% cov.mat.inv
  }

  U <- alpha %*% t(alpha) - cov.mat.inv

  for (i in 1:num.hps) {
    for (j in i:num.hps) {
      L <- sum(cov.mat.inv * (grad_k_alpha[, i] %*% t(grad_k_alpha[, j])))
      M <- sum(t(grad.k.k.inv[, , i]) * grad.k.k.inv[, , j])
      N <- sum(U * cov.mat.hess.array[, , i, j])
      log.marginal.likelihood.hess[i, j] <-
        log.marginal.likelihood.hess[j, i] <- 1/2 * (N + M) - L
    }
  }

  #V <- cov.mat.hess.array
  #
  #for (i in 1:num.hps) {
  #  for (j in 1:num.hps) {
  #    V[, , i, j] <- grad_k_alpha[, i] %*% t(grad_k_alpha[, j])
  #  }
  #}
  # The * rather than %*% in what follows is not a mistake - tr(t(A)%*%B) is
  # the sum of the element-wise product of A and B.
  #temp.multiplicand <- t((alpha %*% t(alpha) - cov.mat.inv))

  #temp.multiplicand <- array(rep(temp.multiplicand, length(hyper.params) + 1),
  #                           dim(cov.mat.grad.array))

  #log.marginal.likelihood.grad <- colSums(1/2 * temp.multiplicand * cov.mat.grad.array, dims=2)

  return(log.marginal.likelihood.hess)
}

get.marginal.likelihood.optimx <- function(par, gp, prior=get.uniform.prior()) {
  -get.marginal.likelihood(gp=gp, sigma.n=par[1], hyper.params=par[-1]) - prior$log.prior(par)
}

get.marginal.likelihood.grad.optimx <- function(par, gp, prior=get.uniform.prior()) {
  -get.marginal.likelihood.grad(gp=gp, sigma.n=par[1], hyper.params=par[-1]) - prior$log.prior.grad(par)
}

optimize.gp.params <- function(gp,
                               init.params,
                               optimx.method,
                               max.iterations,
                               lower,
                               optimx.starttests,
                               optimx.trace,
                               prior=get.uniform.prior()) {
  return(optimx(init.params, get.marginal.likelihood.optimx,
                gr=get.marginal.likelihood.grad.optimx,
                method=optimx.method,
                itnmax=max.iterations,
                lower=lower,
                control=list(starttests=optimx.starttests,
                             trace=optimx.trace,
                             kkt=FALSE),
                gp=gp,
                prior=prior))
}

#' Fit hyperparameters for a GP using maximum likelihood optimisation
#'
#' @param gp a GaussianProcess object
#' @param sigma.n.init an initial value for sigma.n, or NA for random
#' initialisation
#' @param hyper.params.init a vector of initial values for the kernel
#' hyperparameters, or NA for random initialisation. Partial random starts
#' can be specified by setting only some entries in the vector to be NA, or
#' a complete random start can be specified by passing a single NA value.
#' @param random.init.gridsize the number of random initialisation points to
#' test
#' @param resample.top.inits boolean indicating whether to do a second random
#' search based on the (smoothed) distribution of the top results in the first
#' pass
#' @param resample.from number of points to build the distribution from in the
#' second random search
#' @param num.resamples number of samples for the second search
#' @param additional.optimx.runs the number of additional start points to attempt
#' optimisation from
#' @param additional.par.perc.diff the percentage difference required between
#' any two starting parameter vectors
#' @param random.init.scale a scaling factor for the random generation of
#' initialisation points
#' @param abs.min.sigma.n minimum value of sigma.n
#' @param abs.min.hyper.params a vector of minimum values of the hyper params
#' @param optimx.method the optimx method to use for optimisation
#' @param max.iterations max iterations of each optimx run
#' @param optimx.starttests run optimx start tests
#' @param optimx.trace optimx trace level
#' @param verbose verbose output to screen
#'
#' @return a trained GaussianProcess object
#' @export
#'
#' @importFrom cacheMan cached_call
#'
#' @examples
fit.hyperparams <- function(gp,
                            sigma.n.init=NA,
                            hyper.params.init=NA,
                            prior=get.uniform.prior(),
                            random.init.gridsize=50,
                            resample.top.inits=TRUE,
                            resample.from=NULL,
                            num.resamples=NULL,
                            additional.optimx.runs=1,
                            additional.par.perc.diff=0.1,
                            random.init.scale=1,
                            abs.min.sigma.n=0.001,
                            abs.min.hyper.params=0,
                            optimx.method="L-BFGS-B", max.iterations=10000,
                            optimx.starttests=FALSE, optimx.trace=0,
                            verbose=FALSE) {
  if (length(hyper.params.init) == 1 & is.na(hyper.params.init)){
    hyper.params.init <- rep(NA, length(gp$kernel$hyperparam_names))
    names(hyper.params.init) <- gp$kernel$hyperparam_names
  }

  hyper.params.init <- check.hyperparams(gp, hyper.params.init)

  # generate random.init.gridsize random pars and find the best
  best.par <- NULL
  best.log.likelihood <- -Inf

  par.mat <- matrix(0, nrow=random.init.gridsize, ncol=length(hyper.params.init) + 1)
  colnames(par.mat) <- c("sigma.n", names(hyper.params.init))
  log.likelihood.vec <- rep(-Inf, random.init.gridsize)

  successful_init_params <- 0

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
      log.likelihood.temp <- -get.marginal.likelihood.optimx(par=par, gp=gp, prior=prior)

      par.mat[i, ] <- par
      log.likelihood.vec[i] <- log.likelihood.temp

      if (log.likelihood.temp > best.log.likelihood | is.null(best.log.likelihood)) {
        best.log.likelihood <- log.likelihood.temp
        best.par <- par
      }
      successful_init_params <- successful_init_params + 1
    },
    error=function(e) {
      warning(paste("Bad random hyperparams encountered. Error:", e))
    }
    )
  }

  if (successful_init_params == 0 & random.init.gridsize > 0) {
    stop("No initial parameters successfully trained a GP. Check warnings for suppressed error messages")
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
        log.likelihood.temp <- -get.marginal.likelihood.optimx(par=par, gp=gp, prior=prior)

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
  optimx.obj <- optimize.gp.params(gp,
                                   best.par,
                                   optimx.method,
                                   max.iterations,
                                   lower,
                                   optimx.starttests,
                                   optimx.trace,
                                   prior=prior)

  if (optimx.obj["convcode"] > 0) {
    warning("Did not converge!")
    best.log.likelihood <- -Inf
  } else {
    opt.pars <- as.numeric(optimx.obj)[1:(length(hyper.params.init) + 1)]
    names(opt.pars) <- c("sigma.n", names(hyper.params.init))
    best.log.likelihood <- -get.marginal.likelihood.optimx(par=opt.pars, gp=gp, prior=prior)
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

    new.optimx.obj <- optimize.gp.params(gp,
                                         new.par,
                                         optimx.method,
                                         max.iterations,
                                         lower,
                                         optimx.starttests,
                                         optimx.trace,
                                         prior=prior)

    if (new.optimx.obj["convcode"] > 0) {
      warning("Did not converge!")
    } else {
      trained.pars <- as.numeric(new.optimx.obj)[1:(length(hyper.params.init) + 1)]
      names(trained.pars) <- c("sigma.n", names(hyper.params.init))
      new.log.likelihood <- -get.marginal.likelihood.optimx(par=trained.pars, gp=gp, prior=prior)
      if (new.log.likelihood > best.log.likelihood) {
        best.log.likelihood <- new.log.likelihood
        opt.pars <- trained.pars
        optimx.obj <- new.optimx.obj
      }
    }

  }
  if (best.log.likelihood != -Inf) {
    gp$optimized.hyperparams <- opt.pars[-1]
    gp$optimized.sigma.n <- opt.pars[1]
    if (verbose) print(paste("Final.pars after optimisation:"))
    if (verbose) print(opt.pars)
    if (verbose) print(paste("Log likelihood for best starting par optimisation:", best.log.likelihood))
  }

  return(list(gp=gp, optimx.obj=optimx.obj))
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
#' @importFrom cacheMan cached_call
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
  if (!"optimized.hyperparams" %in% names(gp.obj)) {
    stop("No optimized hyperparameters found. This GP has not yet been trained.")
  }

  if (!is.matrix(data)) {
    data <- matrix(data, nrow=length(data))
  }

  num.data.points <- nrow(data)
  kernel.func.list <- list()
  K.star.inst.list <- list()

  K.star <- get_kstar_matrix(kernel=gp.obj$kernel,
                             data.to.predict=data,
                             training.data=gp.obj$training.points,
                             hyper.params=gp.obj$optimized.hyperparams)

  alpha <- cached_call(get.alpha,
                       kernel=gp.obj$kernel,
                       training.points=gp.obj$training.points,
                       training.point.values=gp.obj$training.point.values,
                       sigma.n=gp.obj$optimized.sigma.n,
                       hyper.params=gp.obj$optimized.hyperparams,
                       cache=gp.obj$cache)

  pred.out <- list()
  pred.out[["mean"]] <- K.star %*% alpha

  pred.out[["var"]] <- matrix(NA, nrow=num.data.points, ncol=1)

  cov.mat.chol <- cached_call(get.cov.mat.chol,
                              kernel=gp.obj$kernel,
                              x=gp.obj$training.points,
                              sigma.n=gp.obj$optimized.sigma.n,
                              hyper.params=gp.obj$optimized.hyperparams,
                              cache=gp.obj$cache)

  cov.mat.L.inv <- cached_call(solve,
                               t(cov.mat.chol),
                               cache=gp.obj$cache)

  for (i in seq(num.data.points)) {
    v <- cov.mat.L.inv %*% t(K.star[i, , drop=F])
    pred.out$var[i, 1] <- (get_covariance_matrix(kernel=gp.obj$kernel,
                                                 x=data[i, , drop=FALSE],
                                                 sigma.n=gp.obj$optimized.sigma.n,
                                                 hyper.params=gp.obj$optimized.hyperparams)
                           - t(v) %*% v)
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
#' @param col.index if the predictor data is multi-dimensional, specify which dimension to plot. Other dimensions will be set to their mean.
#' @param ... Additional parameters to pass to \code{plot}
#'
#' @return NULL
#' @export
plot.GaussianProcess <- function(gp.obj, xlim=NULL, num.points=1000, col.index=1, ...) {
  if (is.null(dim(gp.obj$training.points))) {
    x.data <- gp.obj$training.points
  } else {
    x.data <- gp.obj$training.points[, col.index]
  }

  if (is.null(xlim)) {
    x.min <- min(x.data)
    x.max <- max(x.data)

    x.range <- x.max - x.min

    x.min <- x.min - x.range / 4
    x.max <- x.max + x.range / 4
  } else {
    x.min = xlim[1]
    x.max = xlim[2]
  }

  x.values <- seq(from=x.min, to=x.max, length.out=num.points)
  if (is.null(dim(gp.obj$training.points))) {
    preds <- predict(gp.obj, x.values)
  } else {
    input.data <- matrix(colMeans(gp.obj$training.points),
                         nrow=length(x.values),
                         ncol=ncol(gp.obj$training.points),
                         byrow=TRUE)
    input.data[, col.index] <- x.values
    preds <- predict(gp.obj, input.data)
  }

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

    points(x.data, gp.obj$training.point.values)
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
  cov.mat <- get_covariance_matrix(kernel, x, sigma.n, hyper.params)
  y <- rmnorm(n=num.functions, mean=rep(0, nrow(cov.mat)), cov.mat)
  return(y)
}

