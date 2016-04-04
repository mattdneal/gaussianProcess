// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utils.hpp"
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

// [[Rcpp::export]]
NumericVector getCovarianceMatrixHessianArray(NumericMatrix x,
                                              std::string k,
                                              NumericVector sigma_n,
                                              NumericVector hyperParams,
                                              List additionalParams
) {
  int numTrainingPoints = x.nrow();

  int numHyperParams = hyperParams.size() + 1;

  NumericVector covMatHessOut(pow(numTrainingPoints, 2) * pow(numHyperParams, 2));

  int dArr[] = {numTrainingPoints, numTrainingPoints, numHyperParams, numHyperParams};

  std::vector<int> d (dArr, dArr + sizeof(dArr) / sizeof(dArr[0]));

  kernHessPtr kernelHess = selectKernelHess(k);

  NumericMatrix tempHess;
  for (int i=0; i<numTrainingPoints; i++) {
    int arrIndArr[] = {i, i, 0, 0};
    std::vector<int> arrInd (arrIndArr, arrIndArr + sizeof(arrIndArr) / sizeof(arrIndArr[0]));
    int index = arrToVecInd(arrInd, d);
    covMatHessOut(index) = 2;
  }

  if (hyperParams.size() > 0) {
    for (int sample1=0; sample1<numTrainingPoints; sample1++) {
      for (int sample2=0; sample2<=sample1; sample2++) {
        tempHess = kernelHess(x(sample1, _), x(sample2, _), hyperParams, additionalParams);
        for (int i=0; i<hyperParams.size(); i++) {
          for (int j=0; j<hyperParams.size(); j++) {
            int arrIndArr1[] = {sample1, sample2, i + 1, j + 1};
            std::vector<int> arrInd1 (arrIndArr1, arrIndArr1 + sizeof(arrIndArr1) / sizeof(arrIndArr1[0]));
            int index1 = arrToVecInd(arrInd1, d);

            int arrIndArr2[] = {sample2, sample1, i + 1, j + 1};
            std::vector<int> arrInd2 (arrIndArr2, arrIndArr2 + sizeof(arrIndArr2) / sizeof(arrIndArr2[0]));
            int index2 = arrToVecInd(arrInd2, d);

            covMatHessOut(index1) = tempHess(i, j);
            covMatHessOut(index2) = tempHess(i, j);
          }
        }
      }
    }
  }

  NumericVector dim = NumericVector::create(numTrainingPoints, numTrainingPoints, numHyperParams, numHyperParams);

  covMatHessOut.attr("dim") = dim;

  return(covMatHessOut);
}

// A function for calling c++ kernels directly from R
//' Select Built-in C++ Kernels by Name
//'
//' Built in kernels include:
//' \itemize{
//'   \item squaredExponential
//'   \item rationalQuadratic
//'   \item periodic
//'   \item constant
//'   \item generalisedLinear
//'   \item oneDLinear
//'   \item changepoint
//'   \item randomForest
//'   \item neuralNetwork
//'   \item generalisedPolynomial
//'   \item polynomial
//'   \item homogeneousPolynomial
//' }
//'
//' @param kernelName the kernel's name (as a string)
//' @param a the first data point
//' @param b the second data point
//' @param hyperParams the kernel's hyperparameters as a named numeric vector
//' @param additionalParams a list of any additional parameters
//'
//' @export
// [[Rcpp::export]]
NumericVector callKernelByString(std::string kernelName,
                                 NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams) {
  kernPtr kernel = selectKernel(kernelName, FALSE);
  return(kernel(a, b, hyperParams, additionalParams));
}


// [[Rcpp::export]]
NumericVector callKernelGradByString(std::string kernelName,
                                     NumericVector a,
                                     NumericVector b,
                                     NumericVector hyperParams,
                                     List additionalParams) {
  kernPtr kernelGrad = selectKernel(kernelName, TRUE);
  return(kernelGrad(a, b, hyperParams, additionalParams));
}

// [[Rcpp::export]]
NumericVector callKernelHessByString(std::string kernelName,
                                     NumericVector a,
                                     NumericVector b,
                                     NumericVector hyperParams,
                                     List additionalParams) {
  kernHessPtr kernelHess = selectKernelHess(kernelName);
  return(kernelHess(a, b, hyperParams, additionalParams));
}
