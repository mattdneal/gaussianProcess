% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covMatUtils.R
\name{get_covariance_matrix.Kernel}
\alias{get_covariance_matrix.Kernel}
\title{Get the Covariance Matrix for a Kernel}
\usage{
\method{get_covariance_matrix}{Kernel}(kernel, x, sigma.n, hyper.params,
  cache = NULL)
}
\arguments{
\item{kernel}{A Kernel object}

\item{x}{A matrix of data to find the covariance matrix of}

\item{sigma.n}{The standard deviation of the noise to add to the kernel - used to regularise the resulting covariance matrix}

\item{hyper.params}{The hyperparameters of \code{kernel}}

\item{additional.params}{Any additional parameters of \code{kernel} (internal use only)}
}
\value{
A covariance matrix
}
\description{
Get the Covariance Matrix for a Kernel
}

