
#' @importFrom numDeriv grad
check.likelihood.grad <- function(gp, sigma.n, hyper.params)  {
  actual_grad <- get.marginal.likelihood.grad(gp, sigma.n, hyper.params)
  numeric_grad <- grad(function(par) get.marginal.likelihood(gp, par[1], par[-1]), c(sigma.n, hyper.params))
  print(actual_grad)
  print(numeric_grad)
  return(actual_grad - numeric_grad)
}

#' @importFrom numDeriv hessian jacobian
check.likelihood.hessian <- function(gp, sigma.n, hyper.params)  {
  actual_hess <- get.marginal.likelihood.hessian(gp, sigma.n, hyper.params)
  numeric_hess_1 <- hessian(function(par) get.marginal.likelihood(gp, par[1], par[-1]), c(sigma.n, hyper.params))
  numeric_hess_2 <- jacobian(function(par) get.marginal.likelihood.grad(gp, par[1], par[-1]), c(sigma.n, hyper.params))
  print(actual_hess)
  print(numeric_hess_1)
  print(numeric_hess_2)
  return(actual_hess - numeric_hess_1)
}
