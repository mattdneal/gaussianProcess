rmse <- function(prediction, actual) {
  return(sqrt(mean((actual - prediction)^2)))
}

rmse_boot <- function(data, indices, actual) {
  return(rmse(data[indices], actual[indices]))
}