#ifndef __KERNELS__   // if x.h hasn't been included yet...
#define __KERNELS__   //   #define this so the compiler knows it has been included

#include <Rcpp.h>
using namespace Rcpp;


// ****************************************************************************
// * Kernel Hyperparameter Names
// ****************************************************************************

CharacterVector getKernelHyperparamNames(std::string kernelName,
                                         List additionalParams);

// ****************************************************************************
// * Squared Exponential
// ****************************************************************************


const CharacterVector squaredExponentialHPs = CharacterVector::create("l");

NumericVector squaredExponentialKernel(NumericVector a,
                                       NumericVector b,
                                       NumericVector hyperParams,
                                       List additionalParams);


NumericVector squaredExponentialKernelGrad(NumericVector a,
                                           NumericVector b,
                                           NumericVector hyperParams,
                                           List additionalParams);

NumericMatrix squaredExponentialKernelHess(NumericVector a,
                                           NumericVector b,
                                           NumericVector hyperParams,
                                           List additionalParams);


// ****************************************************************************
// * ARD
// ****************************************************************************

NumericVector ARDKernel(NumericVector a,
                        NumericVector b,
                        NumericVector hyperParams,
                        List additionalParams);

NumericVector ARDKernelGrad(NumericVector a,
                            NumericVector b,
                            NumericVector hyperParams,
                            List additionalParams);

NumericMatrix ARDKernelHess(NumericVector a,
                            NumericVector b,
                            NumericVector hyperParams,
                            List additionalParams);


// ****************************************************************************
// * Inverse ARD
// ****************************************************************************

NumericVector inverseARDKernel(NumericVector a,
                               NumericVector b,
                               NumericVector hyperParams,
                               List additionalParams);

NumericVector inverseARDKernelGrad(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams);

NumericMatrix inverseARDKernelHess(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams);

// ****************************************************************************
// * Rational Quadratic
// ****************************************************************************

const CharacterVector rationalQuadraticHPs = CharacterVector::create("l", "alpha");


NumericVector rationalQuadraticKernel(NumericVector a,
                                      NumericVector b,
                                      NumericVector hyperParams,
                                      List additionalParams);


NumericVector rationalQuadraticKernelGrad(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams);

NumericMatrix rationalQuadraticKernelHess(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams);


// ****************************************************************************
// * Periodic
// ****************************************************************************



NumericVector periodicKernel(NumericVector a,
                             NumericVector b,
                             NumericVector hyperParams,
                             List additionalParams);


NumericVector periodicKernelGrad(NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams);



// ****************************************************************************
// * Constant
// ****************************************************************************



NumericVector constantKernel(NumericVector a,
                             NumericVector b,
                             NumericVector hyperParams,
                             List additionalParams);


NumericVector constantKernelGrad(NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams);


// ****************************************************************************
// * One-D Linear
// ****************************************************************************



NumericVector oneDLinearKernel(NumericVector a,
                               NumericVector b,
                               NumericVector hyperParams,
                               List additionalParams);


NumericVector oneDLinearKernelGrad(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams);



// ****************************************************************************
// * Generalised Linear
// ****************************************************************************



NumericVector generalisedLinearKernel(NumericVector a,
                                      NumericVector b,
                                      NumericVector hyperParams,
                                      List additionalParams);


NumericVector generalisedLinearKernelGrad(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams);

// ****************************************************************************
// * Generalised Polynomial
// ****************************************************************************



NumericVector generalisedPolynomialKernel(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams);


NumericVector generalisedPolynomialKernelGrad(NumericVector a,
                                              NumericVector b,
                                              NumericVector hyperParams,
                                              List additionalParams);

// ****************************************************************************
// * Polynomial
// ****************************************************************************



NumericVector polynomialKernel(NumericVector a,
                               NumericVector b,
                               NumericVector hyperParams,
                               List additionalParams);


NumericVector polynomialKernelGrad(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams);
// ****************************************************************************
// * Homogeneous Polynomial
// ****************************************************************************



NumericVector homogeneousPolynomialKernel(NumericVector a,
                               NumericVector b,
                               NumericVector hyperParams,
                               List additionalParams);


NumericVector homogeneousPolynomialKernelGrad(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams);

// ****************************************************************************
// * Changepoint
// ****************************************************************************

double sigmoid(double x,
               double changepoint,
               double transitionRate);


NumericVector changepointKernel(NumericVector a,
                                       NumericVector b,
                                       NumericVector hyperParams,
                                       List additionalParams);


NumericVector changepointKernelGrad(NumericVector a,
                                           NumericVector b,
                                           NumericVector hyperParams,
                                           List additionalParams);


// ****************************************************************************
// * Random Forest
// ****************************************************************************


NumericVector randomForestKernel(NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams);


NumericVector randomForestKernelGrad(NumericVector a,
                                     NumericVector b,
                                     NumericVector hyperParams,
                                     List additionalParams);

// ****************************************************************************
// * Neural Network
// ****************************************************************************


NumericVector neuralNetworkKernel(NumericVector a,
                                  NumericVector b,
                                  NumericVector hyperParams,
                                  List additionalParams);


NumericVector neuralNetworkKernelGrad(NumericVector a,
                                      NumericVector b,
                                      NumericVector hyperParams,
                                      List additionalParams);

// ****************************************************************************
// * Kernel Selector
// ****************************************************************************

// This sets up our pointer for kernel functions so we can pass them around later on
typedef NumericVector (*kernPtr)(NumericVector a,
                       NumericVector b,
                       NumericVector hyperParams,
                       List additionalParams);

// This function selects between kernels. Remember to add new kernels in here.
kernPtr selectKernel(std::string kernelName, bool returnGrad);

// This sets up our pointer for kernel Hessian functions so we can pass them around later on
typedef NumericMatrix (*kernHessPtr)(NumericVector a,
                                     NumericVector b,
                                     NumericVector hyperParams,
                                     List additionalParams);

// This function selects between kernels. Remember to add new kernels in here.
kernHessPtr selectKernelHess(std::string kernelName);

#endif
