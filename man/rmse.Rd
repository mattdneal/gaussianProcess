\name{rmse}
\alias{rmse}
\alias{rmse_boot}
\title{Root Mean Squared Error}
\usage{
rmse(prediction, actual)
rmse_boot(data, indices, actual)
}
\arguments{
  \item{prediction}{ the predicted values. }
  \item{actual}{ the actual/observed values. }
  \item{data}{ the predicted values. }
  \item{indices}{ the indices included in this bootstrap sample. }
}
\description{
Calculates the root mean squared error of predicted values. rmse_boot provides an interface for bootstrap samples of the RMSE using boot.
}
\author{Matthew Neal}
\examples{
prediction <- c(1,5,3,2,2,1,7)
actual <- c(5,3,2,2,2,3,7)
rmse(prediction, actual)

indices <- c(1,4,2)
rmse_boot(prediction, indices, actual)
}
\keyword{rmse}
