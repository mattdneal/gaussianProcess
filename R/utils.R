#' Bootify a Statistic
#'
#' Takes a statistic and returns a wrapper of that statistic suitable for use
#' with \code{\link[boot]{boot}}.
#'
#' @param statistic statistic to bootify
#'
#' @return an R function
bootify <- function(statistic) {
  function(data, indices, actual) {
    return(statistic(data[indices], actual[indices]))
  }
}

#' Regression fit metrics
#'
#' @param prediction predicted values
#' @param actual actual values
#' @param data predicted values
#' @param indices indices included in bootstrap sample
#'
#' @return Statistic value (or vector of values)
#' @export
rmse <- function(prediction, actual) {
  return(sqrt(mean((actual - prediction)^2)))
}

#' @rdname rmse
#' @export
rmse_boot <- bootify(rmse)

#' @rdname rmse
#' @export
coef_of_determination <- function(prediction, actual) {
  return(1 - sum((prediction - actual)^2) / sum((actual - mean(actual))^2))
}

#' @rdname rmse
#' @export
coef_of_determination_boot <- function(data, indices, actual) {
  return(coef_of_determination(data[indices], actual[indices]))
}

#' @rdname rmse
#' @export
mae <- function(prediction, actual) {
  return(mean(abs(prediction - actual)))
}

#' @rdname rmse
#' @export
mae_boot <- function(data, indices, actual){
  return(mae(data[indices], actual[indices]))
}

#' @rdname rmse
#' @export
regression_summary_stats <- function(prediction, actual) {
  out <- c()
  for (fun in c("rmse", "coef_of_determination", "mae")) {
    out <- c(out, get(fun)(prediction, actual))
    names(out)[length(out)] <- fun
  }
  return(out)
}

#' @rdname rmse
#' @export
regression_summary_stats_boot <- bootify(regression_summary_stats)




speedtest.solve <- function(gp.obj) {
  R <-  gp.obj$cov.mat.chol
  y <- gp.obj$training.point.values
  z <- forwardsolve(t(R), y)
  x <- backsolve(R, z)
  return(x)
}

speedtest.inv <- function(gp.obj){
  R <-  gp.obj$cov.mat.chol
  y <- gp.obj$training.point.values
  x.prime <- chol2inv(R) %*% y
  return(x.prime)
}

array.to.vec.coords <- function(x, d) {
  x <- x - 1
  out <- 0
  for (k in seq_along(d)) {
    if (k==1) {
      out <- out + x[k]
    } else {
      out <- out + x[k] * prod(d[1:(k-1)])
    }
  }
  return(out + 1)
}
