#' Scale a kernel
#'
#' @param kernel a kernel object
#' @param kernel_class_name a string naming \code{kernel}
#'
#' @return a kernel object
#' @export
scale_kernel <- function(kernel, kernel_class_name) {
  mt <- create.model.tree()
  mt <- add.kernel(mt, kernel_class_name, kernel)

  const_kern <- create.kernel.object("constant",
                                     hyperparam_names=c("sigma_0"))

  mt <- add.kernel(mt, "constant", const_kern)

  mt <- insert.kernel.instance(mt, 1, kernel_class_name, NULL, NULL)
  mt <- insert.kernel.instance(mt, 1, "constant", "*", NULL)

  scaled_kernel <- create.kernel.object.from.model.tree(mt)

  return(scaled_kernel)
}

#' Scale a list of kernel objects
#'
#' @param kernel_list a list of Kernel objects
#'
#' @return a list of kernel objects
#' @export
scale_kernel_list <- function(kernel_list) {
  output <- list()
  for (kernel_name in names(kernel_list)) {
    kernel <- kernel_list[[kernel_name]]
    output[[kernel_name]] <- scale_kernel(kernel, kernel_name)
  }
  return(output)
}

#' Create a list of polynomial kernels
#'
#' @param degrees a vector containing the degrees to include
#' @param homogeneous boolean indicating whether to add homogeneous polynomial kernels
#' @param poly boolean indicating whether to add non-homogeneous polynomial kernels
#' @param generalised boolean indicating whether to add generalised polynomial kernels
#'
#' @return a list of Kernel objects
#' @export
#'
#' @examples
#' mt <- create.model.tree()
#' mt <- add.polynomial.kernels(mt, degrees=c(2,3))
list_polynomial_kernels <- function(degrees=c(1,2),
                                    homogeneous=TRUE,
                                    poly=TRUE,
                                    generalised=TRUE) {
  output <- list()
  for (degree in degrees) {
    if (homogeneous) {
      kernel <- create.kernel.object(kernel="homogeneousPolynomial",
                                     grad_function=NULL,
                                     hyperparam_names=character(0),
                                     additional_params=list(degree=degree))
      output[[paste("homogPoly", degree, sep="_")]] <- kernel
    }
    if (poly) {
      kernel <- create.kernel.object(kernel="polynomial",
                                     grad_function=NULL,
                                     hyperparam_names=c("c"),
                                     additional_params=list(degree=degree))
      output[[paste("poly", degree, sep="_")]] <- kernel
    }
    if (generalised) {
      kernel <- create.kernel.object(kernel="generalisedPolynomial",
                                     grad_function=NULL,
                                     hyperparam_names=c("c", "l"),
                                     additional_params=list(degree=degree))
      output[[paste("genPoly", degree, sep="_")]] <- kernel
    }
  }
  return(output)
}

#' Return a list containing built in kernels.
#'
#' Not all built in kernels are included - some are data-dependent (RF kernel)
#' and others require additional things to be specified like degree or the
#' dimensions of the data they act on.
#'
#' @return A list of built in kernel objects
#' @export
#'
#' @examples
#' built_in_kernels <- list_built_in_kernels()
#' built_in_kernels
list_built_in_kernels <- function() {
  built_in_kernels <- list()

  built_in_kernels$squaredExponential <- create.kernel.object("squaredExponential",
                                                              hyperparam_names=c("l"))

  built_in_kernels$rationalQuadratic <- create.kernel.object("rationalQuadratic",
                                                             hyperparam_names=c("l", "alpha"))

  built_in_kernels$periodic <- create.kernel.object("periodic",
                                                    hyperparam_names=c("l", "p"))

  built_in_kernels$constant <- create.kernel.object("constant",
                                                    hyperparam_names=c("sigma_0"))

  built_in_kernels$neuralNetwork <- create.kernel.object("neuralNetwork",
                                                         hyperparam_names=c("sigma_0", "sigma"))

  built_in_kernels <- c(built_in_kernels,
                        list_polynomial_kernels(degrees=c(1,2),
                                                homogeneous=FALSE,
                                                poly=FALSE,
                                                generalised=TRUE))

  return(built_in_kernels)
}

#' Return a list of built-in kernels scaled with a constant kernel
#'
#' @return a list of Kernel objects
#' @export
list_built_in_kernels_scaled <- function() {
  built_in_kernels <- list_built_in_kernels()
  kernels_to_scale <- c("squaredExponential",
                        "rationalQuadratic",
                        "periodic",
                        "neuralNetwork")
  scaled_built_in_kernels <- scale_kernel_list(built_in_kernels[kernels_to_scale])
  scaled_built_in_kernels <- c(scaled_built_in_kernels,
                               built_in_kernels[!names(built_in_kernels) %in% kernels_to_scale])
  return(scaled_built_in_kernels)
}
