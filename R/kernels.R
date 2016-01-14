builtin.kernels <- c("squaredExponential",
                     "rationalQuadratic",
                     "periodic",
                     "constant"
                     #"oneDLinear"
                     #"changepoint"
)

# If there are no hyperparams, use numeric() rather than c() so the c++
# code works - it expects a vector not NULL.
builtin.kernel.hyperparam.names <- list("squaredExponential"=c("l"),
                                        "rationalQuadratic"=c("l", "alpha"),
                                        "periodic"=c("l", "p"),
                                        "constant"=c("sigma_0")
                                        #"oneDLinear"=c("intercept","sigma_1")
                                        #"changepoint"=c("changepoint", "transitionRate")
                                        )

builtin.kernel.additional.hyperparams <- list("squaredExponential"=list(),
                                              "rationalQuadratic"=list(),
                                              "periodic"=list(),
                                              "constant"=list()
                                              #"oneDLinear"=list(dimensionIndex=1)
                                              #"changepoint"=list(dimensionIndex=1)
                                             )

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
    print(k(a, b, hyper.params, additional.params))
    print(kg)
    print(kgn)
    if (FALSE & (!isTRUE(all.equal(kg, kgn)) | max.diff > 0.01)) {
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

random.partition.kernel <- function(a, b, hyper.params=NULL, additional.params) {
  partition.function <- additional.params$partitionFunction
  a.partitions <- partition.function(a)
  b.partitions <- partition.function(b)
  num.comparisons <- length(a.partitions)
  return(sum(a.partitions == b.partitions) / num.comparisons)
}

# Random forest kernel
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
