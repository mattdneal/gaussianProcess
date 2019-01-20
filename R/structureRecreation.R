#' Test the ability of the BIC search to identify the data structure
#'
#'
#' @export
test.structure.recreation <- function(base.generator.model,
                                      base.search.model,
                                      num.expansions=1,
                                      hyper.param.range=c(0, 10),
                                      sigma.n.range=c(0.1, 1),
                                      xlim=c(-10, 10),
                                      function.points=1000,
                                      max.points.sampled=50,
                                      min.points.sampled=50,
                                      sample.size.ratio=2,
                                      max.new.nodes=2) {
  model.list <- list()
  next.models.to.expand <- list()

  first.model.name <- as.character(base.generator.model)

  if (is.null(first.model.name)) {
    first.model.name <- "Empty"
  } else {
    model.list[[first.model.name]] <- clone.env(base.generator.model)
  }

  next.models.to.expand[[first.model.name]] <- clone.env(base.generator.model)

  for (i in seq_len(num.expansions)){
    models.to.expand <- next.models.to.expand
    next.models.to.expand <- list()
    for (model in models.to.expand) {
      new.models <- generate.next.models(model)
      model.list <- c(model.list, new.models)
      next.models.to.expand <- c(next.models.to.expand, new.models)
    }
  }

  print(paste(length(model.list), "models generated."))

  generator.model <- model.list[[sample(length(model.list), 1)]]

  for (i in seq_along(generator.model$all.hyper.params)) {
    if (is.na(generator.model$all.hyper.params[i])) {
      generator.model$all.hyper.params[i] <- runif(1, min=hyper.param.range[1], max=hyper.param.range[2])
    }
  }

  sigma.n <- runif(1, min=sigma.n.range[1], max=sigma.n.range[2])
  print("********************************************************************************")
  print("Selected generator model")
  print("********************************************************************************")
  print(generator.model)
  print(generator.model$all.hyper.params)
  print(paste("sigma_n:", sigma.n))

  x <- sort(runif(function.points, min=xlim[1], max=xlim[2]))
  y <- sample.functions.from.model.tree(x, generator.model, sigma.n)

  sample.size <- max.points.sampled
  while(sample.size >= min.points.sampled) {
    plot(x, y, type="l")
    indices <- sample(function.points, sample.size)
    x.sample <- x[indices]
    y.sample <- y[indices]
    points(x.sample, y.sample, col="red")

    gp.searched <- bic.model.search(x.sample, y.sample,
                                    base.search.model,
                                    abs.min.sigma.n=sigma.n.range[1] / 2,
                                    random.init.gridsize=500,
                                    optimx.starttests=FALSE,
                                    additional.optimx.runs=5,
                                    plot.gp=T, verbose=FALSE,
                                    reset.params=TRUE,
                                    max.new.nodes=max.new.nodes)

    y.pred <- predict(gp.searched, x)$mean
    rmse <- sqrt(sum((y - y.pred)^2)/function.points)

    print("********************************************************************************")
    print("Selected generator model")
    print("********************************************************************************")
    print(generator.model)
    print(generator.model$all.hyper.params)
    print(paste("sigma_n:", sigma.n))

    print("********************************************************************************")
    print("Recovered model")
    print("********************************************************************************")
    print(gp.searched$model.tree)
    print(gp.searched$model.tree$all.hyper.params)
    print(paste("sigma_n:", gp.searched$saved.sigma.n))
    print(paste("RMSE:", rmse))
    sample.size <- ceiling(sample.size / sample.size.ratio)
  }
}
