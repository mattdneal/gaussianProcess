% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bootify}
\alias{bootify}
\title{Bootify a Statistic}
\usage{
bootify(statistic)
}
\arguments{
\item{statistic}{statistic to bootify}
}
\value{
an R function
}
\description{
Takes a statistic and returns a wrapper of that statistic suitable for use
with \code{\link[boot]{boot}}.
}

