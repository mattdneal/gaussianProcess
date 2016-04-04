Kernel_class_name <- "Kernel"

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

#' Create a numerical Hessian function for a kernel using the \code{\link[numDeriv]{hessian}} function.
#'
#' @param k kernel function
#' @param method method for creating numeric Hessian, see \code{\link[numDeriv]{hessian}}
#' @param method.args method arguments for creating numeric Hessian, see \code{\link[numDeriv]{hessian}}
#' @param additional.params additional parameters for the kernel function \code{k}
#'
#' @return A function which gives the numeric Hessian of \code{k} w.r.t. the kernel hyperparameters.
#' @export
#'
#' @importFrom numDeriv hessian
create.numeric.hessian <- function(k, method="Richardson", method.args=NULL, additional.params=list()) {
  k.hess <- function(a, b, hyper.params, additional.params) {
    out <- hessian(function(inner.hyper.params) {k(a, b, inner.hyper.params, additional.params)},
                hyper.params,
                method=method,
                method.args=method.args)
    names(out) <- names(hyper.params)
    return(out)
  }
  return(k.hess)
}

#' Create a numerical Hessian function from a kernel grad using the \code{\link[numDeriv]{jacobian}} function.
#'
#' @param k.grad function giving the gradient of the kernel
#' @param method method for creating numeric Jacobian, see \code{\link[numDeriv]{jacobian}}
#' @param method.args method arguments for creating numeric Jacobian, see \code{\link[numDeriv]{jacobian}}
#' @param additional.params additional parameters for the kernel grad function \code{k}
#'
#' @return A function which gives the numeric Hessian for \code{k.grad} w.r.t. the kernel hyperparameters.
#' @export
#'
#' @importFrom numDeriv jacobian
create.numeric.hessian.from.grad <- function(k.grad, method="Richardson", method.args=NULL, additional.params=list()) {
  k.hess <- function(a, b, hyper.params, additional.params) {
    out <- jacobian(function(inner.hyper.params) {k.grad(a, b, inner.hyper.params, additional.params)},
                   hyper.params,
                   method=method,
                   method.args=method.args)
    names(out) <- names(hyper.params)
    return(out)
  }
  return(k.hess)
}

test.kernel.functions <- function(kernel_string,
                                  hyper.param.names,
                                  repetitions=1000,
                                  additional.params=list(),
                                  dimensions=NULL) {
  k <- function(a, b, hyper.params, additional.params) {
    callKernelByString(kernel_string, a, b, hyper.params, additional.params)
  }
  k.grad <- function(a, b, hyper.params, additional.params) {
    callKernelGradByString(kernel_string, a, b, hyper.params, additional.params)
  }
  k.hess <- function(a, b, hyper.params, additional.params) {
    callKernelHessByString(kernel_string, a, b, hyper.params, additional.params)
  }

  k.grad.numeric <- create.numeric.grad(k)
  k.hess.numeric <- create.numeric.hessian.from.grad(k.grad)

  results <- matrix(0, nrow=repetitions, ncol=2)
  colnames(results) <- c("grad.diff", "hess.diff")
  out <- list(params=list(), results=results, outputs=list())

  for (i in 1:repetitions) {
    if (is.null(dimensions)) {
      dims <- ceiling(runif(1) * 100)
    } else {
      dims <- dimensions
    }
    a <- runif(dims) * 100
    b <- runif(dims) * 100
    hyper.params <- rplaw(length(hyper.param.names), -2)
    names(hyper.params) <- hyper.param.names
    kg <- k.grad(a, b, hyper.params, additional.params)
    kgn <- k.grad.numeric(a, b, hyper.params, additional.params)

    kh <- k.hess(a, b, hyper.params, additional.params)
    khn <- k.hess.numeric(a, b, hyper.params, additional.params)

    out$results[i, ] <- c(sum((kg-kgn)^2)^0.5, sum((kh-khn)^2)^0.5)
    out$params[[i]] <- list(a=a, b=b, hyper.params=hyper.params)
    out$outputs[[i]] <- list(k.grad=kg, k.grad.numeric=kgn, k.hess=kh, k.hess.numeric=khn)
  }
  return(out)
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
test.kernel.grad <- function(k,
                             k.grad,
                             hyper.param.names,
                             additional.params,
                             repetitions=1000,
                             dimensions=NULL) {
  k.grad.num <- create.numeric.grad(k, additional.params=additional.params)
  max.diff <- 0
  for (i  in 1:repetitions) {
    if (is.null(dimensions)) {
      dims <- ceiling(runif(1) * 100)
    } else {
      dims <- dimensions
    }
    a <- runif(dims) * 100
    b <- runif(dims) * 100
    hyper.params <- rplaw(length(hyper.param.names), -2)
    names(hyper.params) <- hyper.param.names
    kg <- k.grad(a, b, hyper.params, additional.params)
    kgn <- k.grad.num(a, b, hyper.params, additional.params)
    print(sqrt(sum((kg-kgn)^2))/sqrt(sum((kgn)^2)))
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

#' Create a kernel object
#'
#' Kernel objects provide a unifying wrapper for the various things which
#' can act as kernels (built-in kernels (represented by strings), model trees,
#' and arbitrary R functions).
#'
#' @param kernel_function the kernel function (either a string, a function,
#' or a ModelTree)
#' @param grad_function if \code{kernel_function} is a function, then this
#' should be a function returning the grad w.r.t. the hyperparameters,
#' otherwise this should be null
#' @param hyperparam_names A character vector of the names of the hyperparameters
#' @param additional_params A list of additional parameters
#'
#' @return A Kernel object
#' @export
create.kernel.object <- function(kernel, grad_function=NULL, hess_function=NULL,
                                 hyperparam_names=character(0),
                                 additional_params=list()) {
  kernel_obj <- list()
  kernel_obj$kernel <- kernel

  if (!is.null(grad_function) & !is.function(kernel)) {
    stop("grad_function is not null, but kernel is not a function.")
  }

  if (!is.null(hess_function) & !is.function(kernel)) {
    stop("hess_function is not null, but kernel is not a function.")
  }

  kernel_obj$grad <- NULL
  kernel_obj$hess <- NULL

  if (is.function(kernel)) {
    # Gradient
    if (is.function(grad_function)) {
      # Grad supplied.
      kernel_obj$grad <- grad_function
    } else if (is.null(grad_function)) {
      # No grad supplied. Make a grad.
      kernel_obj$grad <- create.numeric.grad(kernel)
    } else {
      stop("Invalid grad_function")
    }

    # Hessian
    if (is.function(hess_function)) {
      # Hessian supplied.
      kernel_obj$hess <- hess_function
    } else if (is.null(hess_function)) {
      # No Hessian supplied. Make a Hessian function
      if (is.function(grad_function)) {
        kernel_obj$hess <- create.numeric.hessian.from.grad(grad_function)
      } else {
        kernel_obj$hess <- create.numeric.hessian(kernel)
      }
    } else {
      stop("Invalid hess_function")
    }
  }

  if (is.character(kernel)) {
    hyperparam_names <- getKernelHyperparamNames(kernel, additional_params)
  }

  if (length(hyperparam_names) == 0) {
    hyperparam_names <- character(0)
  }

  kernel_obj$hyperparam_names <- hyperparam_names

  kernel_obj$additional_params <- additional_params
  class(kernel_obj) <- Kernel_class_name

  return(kernel_obj)
}

#' Print a kernel
#'
#' @param kernel
#'
#' @return NULL
#' @export
print.Kernel <- function(kernel) {
  if (is.character(kernel$kernel)) {
    cat(paste('  - Built-in kernel: ', kernel$kernel, '\n', sep=""))
  }
  if (length(kernel$hyperparam_names) > 0) {
    cat(paste('  - Hyperparameters:\n', sep=""))
    for (param in kernel$hyperparam_names) {
      cat(paste('    - ', param, '\n', sep=""))
    }
  } else {
    cat(paste('  - No hyperparameters\n', sep=""))
  }
  if (length(kernel$additional_params) > 0) {
    cat(paste('  - Additional Parameters:\n', sep=""))
    for (addparam in names(kernel$additional_params)) {
      cat(paste('    - ', addparam, '\n', sep=""))
    }
  }else {
    cat(paste('  - No additional parameters \n', sep=""))
  }
}

#' Create Kernel object from model tree
#'
#' @param model_tree
#'
#' @return a kernel object
#' @export
create.kernel.object.from.model.tree <- function(model_tree) {
  hyperparam_names <- names(model_tree$all.hyper.params)
  kernel <- create.kernel.object(kernel=model_tree,
                                 grad_function=NULL,
                                 hyperparam_names=hyperparam_names,
                                 additional_params=list()
                                 )
  return(kernel)
}

#' Create an ARD kernel object
#'
#' @param dimensions
#'
#' @return a kernel object
#' @export
create.ard.kernel <- function(dimensions, inverse=T) {
  if (inverse) {
    k <- "inverseARD"
  } else {
    k <- "ARD"
  }
  hyperparam_names <- paste("l", 1:dimensions, sep="")
  kernel <- create.kernel.object(k, grad_function=NULL, hyperparam_names=hyperparam_names)
}

create.gen.nn.kernel <- function(dimensions) {
  k <- "generalisedNeuralNetwork"
  hyperparam_names <- paste("sigma", 0:(dimensions + 1), sep="")
  kernel <- create.kernel.object(k, grad_function=NULL, hyperparam_names=hyperparam_names)
}

set.kernel.hash <- function(kernel, cache) {
  attr(kernel, cacheMan:::hash_attr) <- cacheMan::hash_object(kernel, cache)
  if (class(kernel) == Kernel_class_name) {
    kernel$kernel <- set.kernel.hash(kernel$kernel, cache)
  }

  if (class(kernel) == ModelTree_class_name) {
    for (i in seq_along(kernel$kernel.objects)) {
      kernel$kernel.objects[[i]] <- set.kernel.hash(kernel$kernel.objects[[i]], cache)
    }
  }
  return(kernel)
}
