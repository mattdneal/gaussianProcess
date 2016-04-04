#ifndef __COVMATUTILS__   // if x.h hasn't been included yet...
#define __COVMATUTILS__   //   #define this so the compiler knows it has been included

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

/* Because of the cost of context switching to call the R function k, this is actually
 slower than doing it all in R.  */
NumericMatrix getCovarianceMatrixCpp(NumericMatrix x,
                                     Function k,
                                     NumericVector sigma_n,
                                     NumericVector hyperParams);

NumericMatrix getCovarianceMatrixBuiltInCpp(NumericMatrix x,
                                            std::string k,
                                            NumericVector sigma_n,
                                            NumericVector hyperParams,
                                            List additionalParams);

NumericVector getCovarianceMatrixGradArray(NumericMatrix x,
                                           std::string k,
                                           NumericVector sigma_n,
                                           NumericVector hyperParams,
                                           List additionalParams
                                           );


// A function for calling c++ kernels directly from R
NumericVector callKernelByString(std::string kernelName,
                                 NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams);


// A function for calling c++ kernels directly from R
NumericVector callKernelGradByString(std::string kernelName,
                                     NumericVector a,
                                     NumericVector b,
                                     NumericVector hyperParams,
                                     List additionalParams);

// A function for calling c++ kernels directly from R
NumericVector callKernelHessByString(std::string kernelName,
                                     NumericVector a,
                                     NumericVector b,
                                     NumericVector hyperParams,
                                     List additionalParams);

#endif
