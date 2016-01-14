// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "covMatUtils.hpp"
#include "kernels.hpp"
using namespace Rcpp;

// Because of the cost of context switching to call the R function k, this is actually
// slower than doing it all in R. 
// [[Rcpp::export]]
NumericMatrix getCovarianceMatrixCpp(NumericMatrix x, 
                                     Function k, 
                                     NumericVector sigma_n, 
                                     NumericVector hyperParams) {
  int numTrainingPoints = x.nrow();
  
  NumericMatrix covMat(numTrainingPoints, numTrainingPoints);
  
  
  for (int sample1=0; sample1<numTrainingPoints; sample1++) {
    for (int sample2=0; sample2<sample1; sample2++) {
      // Feels like there should be a better way of doing this
      NumericVector covariance(k(x(sample1, _), x(sample2, _), hyperParams));
      covMat(sample1, sample2) = covariance(0);
      covMat(sample2, sample1) = covariance(0);
    }
  }
  for (int i=0; i<numTrainingPoints; i++) {
    covMat(i, i) = covMat(i, i) + pow(sigma_n(0), 2);
  }
  return(covMat);
}


// [[Rcpp::export]]
NumericMatrix getCovarianceMatrixBuiltInCpp(NumericMatrix x, 
                                            std::string k, 
                                            NumericVector sigma_n, 
                                            NumericVector hyperParams,
                                            List additionalParams) {
  int numTrainingPoints = x.nrow();
  
  NumericMatrix covMat(numTrainingPoints, numTrainingPoints);
  
  kernPtr kernel = selectKernel(k, FALSE);
  
  double covariance;
  
  for (int sample1=0; sample1<numTrainingPoints; sample1++) {
    for (int sample2=0; sample2<=sample1; sample2++) {
      covariance = kernel(x(sample1, _), x(sample2, _), hyperParams, additionalParams)(0);
      covMat(sample1, sample2) = covariance;
      covMat(sample2, sample1) = covariance;
    }
  }
  
  for (int i=0; i<numTrainingPoints; i++) {
    covMat(i, i) = covMat(i, i) + pow(sigma_n(0), 2);
  }
  return(covMat);
}

  
// [[Rcpp::export]]
NumericVector getCovarianceMatrixGradArray(NumericMatrix x,
                                           std::string k,
                                           NumericVector sigma_n, 
                                           NumericVector hyperParams,
                                           List additionalParams
                                           ) {
  int numTrainingPoints = x.nrow();
  
  arma::cube covMatGrad(numTrainingPoints, 
                        numTrainingPoints, 
                        hyperParams.size() + 1,
                        arma::fill::zeros
                          );
  
  kernPtr kernelGrad = selectKernel(k, TRUE);
  
  NumericVector tempGrad;
  for (int i=0; i<numTrainingPoints; i++) {
    covMatGrad(i, i, 0) = 2 * sigma_n(0);
  }  
  
  if (hyperParams.size() > 0) {
    for (int sample1=0; sample1<numTrainingPoints; sample1++) {
      for (int sample2=0; sample2<=sample1; sample2++) {
        tempGrad = kernelGrad(x(sample1, _), x(sample2, _), hyperParams, additionalParams);
        for (int i=0; i<hyperParams.size(); i++) {
          covMatGrad(sample1, sample2, i+1) = tempGrad(i); 
          covMatGrad(sample2, sample1, i+1) = tempGrad(i); 
        }
      }
    }
  }
  return(wrap(covMatGrad));
}


// A function for calling c++ kernels directly from R
// [[Rcpp::export]]
NumericVector callKernelByString(std::string kernelName, 
                                 NumericVector a, 
                                 NumericVector b, 
                                 NumericVector hyperParams, 
                                 List additionalParams) {
  kernPtr kernel = selectKernel(kernelName, FALSE);
  return(kernel(a, b, hyperParams, additionalParams));
}

/*** R
*/
