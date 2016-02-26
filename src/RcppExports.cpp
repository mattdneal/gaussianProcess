// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// getCovarianceMatrixCpp
NumericMatrix getCovarianceMatrixCpp(NumericMatrix x, Function k, NumericVector sigma_n, NumericVector hyperParams);
RcppExport SEXP gaussianProcess_getCovarianceMatrixCpp(SEXP xSEXP, SEXP kSEXP, SEXP sigma_nSEXP, SEXP hyperParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< Function >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma_n(sigma_nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    __result = Rcpp::wrap(getCovarianceMatrixCpp(x, k, sigma_n, hyperParams));
    return __result;
END_RCPP
}
// getCovarianceMatrixBuiltInCpp
NumericMatrix getCovarianceMatrixBuiltInCpp(NumericMatrix x, std::string k, NumericVector sigma_n, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_getCovarianceMatrixBuiltInCpp(SEXP xSEXP, SEXP kSEXP, SEXP sigma_nSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma_n(sigma_nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(getCovarianceMatrixBuiltInCpp(x, k, sigma_n, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// getCovarianceMatrixGradArray
NumericVector getCovarianceMatrixGradArray(NumericMatrix x, std::string k, NumericVector sigma_n, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_getCovarianceMatrixGradArray(SEXP xSEXP, SEXP kSEXP, SEXP sigma_nSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma_n(sigma_nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(getCovarianceMatrixGradArray(x, k, sigma_n, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// callKernelByString
NumericVector callKernelByString(std::string kernelName, NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_callKernelByString(SEXP kernelNameSEXP, SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::string >::type kernelName(kernelNameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(callKernelByString(kernelName, a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// squaredExponentialKernel
NumericVector squaredExponentialKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_squaredExponentialKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(squaredExponentialKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// squaredExponentialKernelGrad
NumericVector squaredExponentialKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_squaredExponentialKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(squaredExponentialKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// ARDKernel
NumericVector ARDKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_ARDKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(ARDKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// ARDKernelGrad
NumericVector ARDKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_ARDKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(ARDKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// inverseARDKernel
NumericVector inverseARDKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_inverseARDKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(inverseARDKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// inverseARDKernelGrad
NumericVector inverseARDKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_inverseARDKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(inverseARDKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// rationalQuadraticKernel
NumericVector rationalQuadraticKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_rationalQuadraticKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(rationalQuadraticKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// rationalQuadraticKernelGrad
NumericVector rationalQuadraticKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_rationalQuadraticKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(rationalQuadraticKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// periodicKernel
NumericVector periodicKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_periodicKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(periodicKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// periodicKernelGrad
NumericVector periodicKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_periodicKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(periodicKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// constantKernel
NumericVector constantKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_constantKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(constantKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// constantKernelGrad
NumericVector constantKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_constantKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(constantKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// oneDLinearKernel
NumericVector oneDLinearKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_oneDLinearKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(oneDLinearKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// oneDLinearKernelGrad
NumericVector oneDLinearKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_oneDLinearKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(oneDLinearKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// generalisedLinearKernel
NumericVector generalisedLinearKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_generalisedLinearKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(generalisedLinearKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// generalisedLinearKernelGrad
NumericVector generalisedLinearKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_generalisedLinearKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(generalisedLinearKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// generalisedPolynomialKernel
NumericVector generalisedPolynomialKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_generalisedPolynomialKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(generalisedPolynomialKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// generalisedPolynomialKernelGrad
NumericVector generalisedPolynomialKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_generalisedPolynomialKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(generalisedPolynomialKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// polynomialKernel
NumericVector polynomialKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_polynomialKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(polynomialKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// polynomialKernelGrad
NumericVector polynomialKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_polynomialKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(polynomialKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// homogeneousPolynomialKernel
NumericVector homogeneousPolynomialKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_homogeneousPolynomialKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(homogeneousPolynomialKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// homogeneousPolynomialKernelGrad
NumericVector homogeneousPolynomialKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_homogeneousPolynomialKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(homogeneousPolynomialKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// changepointKernel
NumericVector changepointKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_changepointKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(changepointKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// changepointKernelGrad
NumericVector changepointKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_changepointKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(changepointKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// randomForestKernel
NumericVector randomForestKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_randomForestKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(randomForestKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// randomForestKernelGrad
NumericVector randomForestKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_randomForestKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(randomForestKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// neuralNetworkKernel
NumericVector neuralNetworkKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_neuralNetworkKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(neuralNetworkKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// neuralNetworkKernelGrad
NumericVector neuralNetworkKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_neuralNetworkKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(neuralNetworkKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// generalNeuralNetworkKernel
NumericVector generalNeuralNetworkKernel(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_generalNeuralNetworkKernel(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(generalNeuralNetworkKernel(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// generalNeuralNetworkKernelGrad
NumericVector generalNeuralNetworkKernelGrad(NumericVector a, NumericVector b, NumericVector hyperParams, List additionalParams);
RcppExport SEXP gaussianProcess_generalNeuralNetworkKernelGrad(SEXP aSEXP, SEXP bSEXP, SEXP hyperParamsSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type hyperParams(hyperParamsSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(generalNeuralNetworkKernelGrad(a, b, hyperParams, additionalParams));
    return __result;
END_RCPP
}
// getTreeHeight
NumericVector getTreeHeight(NumericMatrix tree);
RcppExport SEXP gaussianProcess_getTreeHeight(SEXP treeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type tree(treeSEXP);
    __result = Rcpp::wrap(getTreeHeight(tree));
    return __result;
END_RCPP
}
// navigateRFTree
NumericVector navigateRFTree(NumericMatrix tree, NumericVector data, IntegerVector isFactor, NumericVector inputHeight);
RcppExport SEXP gaussianProcess_navigateRFTree(SEXP treeSEXP, SEXP dataSEXP, SEXP isFactorSEXP, SEXP inputHeightSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type isFactor(isFactorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type inputHeight(inputHeightSEXP);
    __result = Rcpp::wrap(navigateRFTree(tree, data, isFactor, inputHeight));
    return __result;
END_RCPP
}
// kahanSum
double kahanSum(NumericVector summands);
RcppExport SEXP gaussianProcess_kahanSum(SEXP summandsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type summands(summandsSEXP);
    __result = Rcpp::wrap(kahanSum(summands));
    return __result;
END_RCPP
}
// sumSQuaredDiffs
double sumSQuaredDiffs(NumericVector a, NumericVector b);
RcppExport SEXP gaussianProcess_sumSQuaredDiffs(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    __result = Rcpp::wrap(sumSQuaredDiffs(a, b));
    return __result;
END_RCPP
}
// sumSQuaredDiffsPartial
double sumSQuaredDiffsPartial(NumericVector a, NumericVector b, List additionalParams);
RcppExport SEXP gaussianProcess_sumSQuaredDiffsPartial(SEXP aSEXP, SEXP bSEXP, SEXP additionalParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< List >::type additionalParams(additionalParamsSEXP);
    __result = Rcpp::wrap(sumSQuaredDiffsPartial(a, b, additionalParams));
    return __result;
END_RCPP
}
// bitWiseAnd
int bitWiseAnd(NumericVector a, NumericVector b);
RcppExport SEXP gaussianProcess_bitWiseAnd(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    __result = Rcpp::wrap(bitWiseAnd(a, b));
    return __result;
END_RCPP
}
