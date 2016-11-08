# Priors

prior_class_name <- "HyperParamPrior"

#' Return an object representing a uniform prior over hyperparameters
#'
#' @return a HyperParamPrior object
#' @export
#'
#' @examples
get.uniform.prior <- function() {

  uniform.prior <- list()
  class(uniform.prior) <- prior_class_name

  uniform.prior$log.prior <- function(hyperparams) {
    return(0)
  }

  uniform.prior$log.prior.grad <- function(hyperparams) {
    return(rep(0, length(hyperparams)))
  }

  uniform.prior$log.prior.hessian <- function(hyperparams) {
    hessian <- matrix(0, length(hyperparams), length(hyperparams))
    colnames(hessian) <- names(hyperparams)
    rownames(hessian) <- names(hyperparams)
    return(hessian)
  }

  return(uniform.prior)
}

#' Return an object representing a normal prior over hyperparameters
#'
#' @param sigma the standard deviation of the prior
#'
#' @return a HyperParamPrior object
#' @export
#'
#' @examples
get.normal.prior <- function(sigma) {

  normal.prior <- list()
  class(normal.prior) <- prior_class_name

  sigma <- abs(sigma)

  normal.prior$sigma <- sigma

  normal.prior$log.prior <- function(hyperparams) {
    return(-length(hyperparams) / 2 * log(sigma * 2 * pi) - sum(hyperparams^2) / (2 * sigma^2))
  }

  normal.prior$log.prior.grad <- function(hyperparams) {
    return(-hyperparams / sigma^2)
  }

  normal.prior$log.prior.hessian <- function(hyperparams) {
    hessian <- diag(-1/sigma^2, length(hyperparams))
    colnames(hessian) <- names(hyperparams)
    rownames(hessian) <- names(hyperparams)
    return(hessian)
  }

  return(normal.prior)
}
