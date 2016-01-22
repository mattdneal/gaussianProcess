
#' Create a numerical gradient function for a kernel using the \code{\link[numDeriv]{grad}} function.
#'
#' @param k kernel function
#' @param method method for creating numeric grad, see \code{\link[numDeriv]{grad}}
#' @param method.args method arguments for creating numeric grad, see \code{\link[numDeriv]{grad}}
#' @param additional.params additional parameters for the kernel function \code{k}
#'
#' @return A function which gives the numeric grad of \code{k} w.r.t. the kernel hyperparameters.
#' @export
#'
#' @importFrom numDeriv grad
create.numeric.grad <- function(k, method="Richardson", method.args=NULL, additional.params=list()) {
  k.grad <- function(a, b, hyper.params, additional.params) {
    out <- grad(function(inner.hyper.params) {k(a, b, inner.hyper.params, additional.params)},
                hyper.params,
                method=method,
                method.args=method.args)
    names(out) <- names(hyper.params)
    return(out)
  }
  return(k.grad)
}

#' Test an Analytic Gradient Against a Numeric Gradient
#'
#' Creates a numerical approximation to the gradient function of a kernel (using \code{\link{create.numeric.grad}})
#' and compares it to an analytic gradient function. If the maximum relative difference exceeds 0.01, returns the failing inputs,
#' otherwise returns TRUE.
#'
#' @param k a kernel function
#' @param k.grad the analytic gradient to test
#' @param hyper.param.names names of the kernel's hyperparameters
#' @param additional.params any additional parameters of the kernel
#' @param repetitions number of repetitions to attempt
#'
#' @return If the gradients match, returns TRUE. If a mismatch is found, returns a list with named entries:
#' \itemize{
#'   \item a, b, hyper.params - the parameters passed to the gradient functions
#'   \item grad - return value of the analytic gradient
#'   \item grad.num - return value of the numerical gradient
#'   \item max.diff - the maximum relative difference observed between the two gradient functions so far (including the failed set of inputs)
#'   \item i - the number of inputs tested (including the failed set of inputs)
#' }
#' @export
test.kernel.grad <- function(k, k.grad, hyper.param.names, additional.params, repetitions=1000) {
  k.grad.num <- create.numeric.grad(k, additional.params=additional.params)
  max.diff <- 0
  for (i  in 1:repetitions) {
    dims <- ceiling(runif(1) * 100)
    a <- runif(dims) * 100
    b <- runif(dims) * 100
    hyper.params <- rplaw(length(hyper.param.names), -2)
    names(hyper.params) <- hyper.param.names
    kg <- k.grad(a, b, hyper.params, additional.params)
    kgn <- k.grad.num(a, b, hyper.params, additional.params)
    #print(sqrt(sum((kg-kgn)^2))/sqrt(sum((kgn)^2)))
    #print(sqrt(sum((a-b)^2)))
    max.diff <- max(max.diff, sqrt(sum((kg-kgn)^2))/sqrt(sum((kgn)^2)))
    #print(kg / kgn)
    #print(k(a, b, hyper.params, additional.params))
    #print(kg)
    #print(kgn)
    if ((max.diff > 0.01)) {
      return(list(a=a,
                  b=b,
                  hyper.params=hyper.params,
                  grad=kg,
                  grad.num=kgn,
                  max.diff=max.diff,
                  i=i
      )
      )
    }
  }
  print(paste("Maximum relative difference:", max.diff))
  return(TRUE)
}

squared.exponential.kernel.ard <- function(a, b, hyper.params) {
  l <- hyper.params[-1]
  if (length(l) != length(a)) {
    stop("Not enough length hyperparameters.")
  }
  return(as.numeric(hyper.params["sigma"]^2 * exp(-1/2 * sum(((a - b) / l)^2))))
}

squared.exponential.kernel.ard.grad <- function(a, b, hyper.params) {
  print("grad called")
  l <- hyper.params[-1]
  if (length(l) != length(a)) {
    stop("Not enough length hyperparameters.")
  }

  precalc.exp <- as.numeric(hyper.params["sigma"] * exp(-1/2 * sum(((a - b) / l)^2)))

  ret <- c(sigma=as.numeric(2 * precalc.exp))

  #for (i in 2:length(hyper.params)) {
  #  ret <- c(ret,
  #           as.numeric(((a[i-1] - b[i-1])^2)/(hyper.params[i]^3) *
  #                        hyper.params["sigma"] *
  #                        precalc.exp)
  #           )
  #  if (i%%1000==0) print(i)
  #}

  ret <- c(ret,
           as.numeric(((a - b)^2)/(l^3) *
                        hyper.params["sigma"] *
                        precalc.exp
           )
  )

  names(ret) <- names(hyper.params)

  return(ret)
}

#' Random Partition Kernel
#'
#' Takes a partition function and returns the covariance between two points
#' derived from that partition function. The partition function should return
#' a vector of integers, with each integer representing the result of a different
#' classifier.
#'
#' @param a first data point
#' @param b second data point
#' @param hyper.params an empy vector (included for compatibility with kernel function format)
#' @param additional.params A list with a named element "partitionFunction" containing the partition function which defines the kernel
#'
#' @return The covariance between \code{a} and \code{b}
#' @export
random.partition.kernel <- function(a, b, hyper.params=NULL, additional.params) {
  partition.function <- additional.params$partitionFunction
  a.partitions <- partition.function(a)
  b.partitions <- partition.function(b)
  num.comparisons <- length(a.partitions)
  return(sum(a.partitions == b.partitions) / num.comparisons)
}

#' Random forest kernel
#'
#' @param rf a randomForest object
#'
#' @return a function
#' @export
#'
#' @importFrom randomForest getTree
#' @importFrom hash hash
#'
#' @seealso \code{\link[randomForest]{randomForest}}
random.forest.partition.function.generator <- function(rf) {
  num.trees <- rf$ntree
  trees <- list()
  heights <- numeric(num.trees)
  cache <- hash()
  for (i in 1:num.trees) {
    trees[[i]] <- getTree(rf, i)
    heights[i] <- sample(0:get.tree.height(trees[[i]]), 1)
  }
  function(data) {
    # Using toString as our hash seems fraught with peril. Discuss.
    key <- toString(data)
    if (has.key(key, cache)) {
      return(cache[[key]])
    } else {
      ret <- sapply(1:num.trees, function(k) navigate.rf.tree(trees[[k]], data, heights[k]))
      cache[[key]] <- ret
      cache <<- cache
      return(ret)
    }
  }
}
