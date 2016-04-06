#define _USE_MATH_DEFINES
#include "utils.hpp"
#include "kernels.hpp"
#include <Rcpp.h>
#include "rfTreeUtils.hpp"
#include <math.h>
#include <string>
using namespace Rcpp;

// ****************************************************************************
// * Kernel Hyperparameter Names
// ****************************************************************************

// [[Rcpp::export]]
CharacterVector getKernelHyperparamNames(std::string kernelName,
                                         List additionalParams) {
  CharacterVector dummy;

  if (kernelName == "squaredExponential") {
    return(squaredExponentialHPs);
  } else if (kernelName == "ARD") {
    return(dummy);
  } else if (kernelName == "inverseARD") {
    return(dummy);
  } else if (kernelName == "rationalQuadratic") {
    return(rationalQuadraticHPs);
  } else if (kernelName == "periodic") {
    return(dummy);
  } else if (kernelName == "constant") {
    return(dummy);
  } else if (kernelName == "generalisedLinear") {
    return(dummy);
  } else if (kernelName == "oneDLinear") {
    return(dummy);
  } else if (kernelName == "changepoint") {
    return(dummy);
  } else if (kernelName == "randomForest") {
    return(dummy);
  } else if (kernelName == "neuralNetwork") {
    return(dummy);
  } else if (kernelName == "generalisedNeuralNetwork") {
    return(dummy);
  } else if (kernelName == "generalisedPolynomial") {
    return(dummy);
  } else if (kernelName == "polynomial") {
    return(dummy);
  } else if (kernelName == "homogeneousPolynomial") {
    return(dummy);
  } else {
    throw std::range_error("Incorrect kernel specified");
  }
}

// ****************************************************************************
// * Squared Exponential
// ****************************************************************************

// [[Rcpp::export]]
NumericVector squaredExponentialKernel(NumericVector a,
                                       NumericVector b,
                                       NumericVector hyperParams,
                                       List additionalParams) {

  checkHPNames(squaredExponentialHPs, hyperParams.names());

  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);

  double l = hyperParams[0];

  double result = exp(-sum / ( 2 * pow(l, 2)));

  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector squaredExponentialKernelGrad(NumericVector a,
                                           NumericVector b,
                                           NumericVector hyperParams,
                                           List additionalParams) {

  checkHPNames(squaredExponentialHPs, hyperParams.names());

  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);

  double l = hyperParams[0];

  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);
  result[0] = sum / pow(l, 3) *
    exp(-sum/(2 * pow(l, 2)));
  return result;
}

// [[Rcpp::export]]
NumericMatrix squaredExponentialKernelHess(NumericVector a,
                                           NumericVector b,
                                           NumericVector hyperParams,
                                           List additionalParams) {

  checkHPNames(squaredExponentialHPs, hyperParams.names());

  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);

  double l = hyperParams[0];

  NumericMatrix result(hyperParams.size());
  result.attr("dimnames") = List::create(squaredExponentialHPs, squaredExponentialHPs);

  result(0, 0) = sum / pow(l, 4) * exp(-sum / ( 2 * pow(l, 2))) * (sum / pow(l, 2) - 3);
  return result;
}


// ****************************************************************************
// * ARD
// ****************************************************************************

// [[Rcpp::export]]
NumericVector ARDKernel(NumericVector a,
                        NumericVector b,
                        NumericVector hyperParams,
                        List additionalParams) {
  double sum = sumSQuaredDiffs(a / hyperParams, b / hyperParams);
  double result = exp(-sum / 2);
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector ARDKernelGrad(NumericVector a,
                            NumericVector b,
                            NumericVector hyperParams,
                            List additionalParams) {
  double sum = sumSQuaredDiffs(a / hyperParams, b / hyperParams);
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);
  for (int i=0; i < result.length(); i++) {
    result[i] = pow(a[i] - b[i], 2) / pow(hyperParams[i], 3) *
      exp(-sum / 2);
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix ARDKernelHess(NumericVector a,
                            NumericVector b,
                            NumericVector hyperParams,
                            List additionalParams) {
  double sum = sumSQuaredDiffs(a / hyperParams, b / hyperParams);
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericMatrix result(hyperParams.size());
  result.attr("dimnames") = List::create(hyperParams.names(), hyperParams.names());

  for (int i=0; i < hyperParams.length(); i++) {
    for (int j=i; j < hyperParams.length(); j++) {
      result(i, j) = pow(a[i] - b[i], 2) / pow(hyperParams[i], 3) *
                      pow(a[j] - b[j], 2) / pow(hyperParams[j], 3) *
                      exp(-sum / 2);
      if (i == j) {
        result(i, j) -= 3 * pow(a[i] - b[i], 2) / pow(hyperParams[i], 4) * exp(-sum / 2);
      }
      result(j, i) = result(i, j);
    }
  }
  return result;
}


// ****************************************************************************
// * Inverse ARD
// ****************************************************************************

// [[Rcpp::export]]
NumericVector inverseARDKernel(NumericVector a,
                        NumericVector b,
                        NumericVector hyperParams,
                        List additionalParams) {
  double sum = sumSQuaredDiffs(a * hyperParams, b * hyperParams);
  double result = exp(-sum / 2);
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector inverseARDKernelGrad(NumericVector a,
                            NumericVector b,
                            NumericVector hyperParams,
                            List additionalParams) {
  double sum = sumSQuaredDiffs(a * hyperParams, b * hyperParams);
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  for (int i=0; i < result.length(); i++) {
    result[i] = -hyperParams[i] * pow(a[i] - b[i], 2) * exp(-sum / 2);
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix inverseARDKernelHess(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams) {
  double sum = sumSQuaredDiffs(a * hyperParams, b * hyperParams);
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericMatrix result(hyperParams.size());
  result.attr("dimnames") = List::create(hyperParams.names(), hyperParams.names());

  for (int i=0; i < hyperParams.length(); i++) {
    for (int j=i; j < hyperParams.length(); j++) {
      result(i, j) = hyperParams[i] * pow(a[i] - b[i], 2) *
        hyperParams[j] * pow(a[j] - b[j], 2) *
        exp(-sum / 2);
      if (i == j) {
        result(i, j) -= exp(-sum / 2);
      }
      result(j, i) = result(i, j);
    }
  }
  return result;
}

// ****************************************************************************
// * Rational Quadratic
// ****************************************************************************

// [[Rcpp::export]]
NumericVector rationalQuadraticKernel(NumericVector a,
                                      NumericVector b,
                                      NumericVector hyperParams,
                                      List additionalParams) {

  checkHPNames(rationalQuadraticHPs, hyperParams.names());

  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);
  double l = hyperParams[0];
  double alpha = hyperParams[1];
  double result = pow(1 + sum / ( 2 * alpha * pow(l, 2)), -alpha);
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector rationalQuadraticKernelGrad(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams) {

  checkHPNames(rationalQuadraticHPs, hyperParams.names());

  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  double l = hyperParams[0];
  double alpha = hyperParams[1];
  //Rcout << sigma << "; " << l << "; " << alpha << ";\n";

  double u = sum / (2 * pow(l, 2));
  double v = 1 + u / alpha;
  double w = sum / pow(l, 3);
  double x = u / (alpha + u);
  //Rcout << u << "; " << v << "; " << w << "; " << x << ";\n";

  result[0] = w * pow(v, -(alpha + 1));
  result[1] = pow(v, -alpha) * (x - log(v));

  return result;
}

// [[Rcpp::export]]
NumericMatrix rationalQuadraticKernelHess(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams) {

  checkHPNames(rationalQuadraticHPs, hyperParams.names());

  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);

  NumericMatrix result(hyperParams.size());
  result.attr("dimnames") = List::create(hyperParams.names(), hyperParams.names());

  double l = hyperParams[0];
  double alpha = hyperParams[1];

  double k = rationalQuadraticKernel(a, b, hyperParams, additionalParams)[0];

  NumericVector kGrad = rationalQuadraticKernelGrad(a, b, hyperParams, additionalParams);

  result(0, 0) = sum / pow(l, 4) *
    pow(1 + sum / (2 * alpha * pow(l, 2)), -(alpha + 1)) *
    (sum *
      (alpha + 1) /
      (alpha * pow(l, 2)) *
      pow(1 + sum / (2 * alpha * pow(l, 2)), -1) -
      3);

  result(0, 1) = kGrad[0] * (sum / (2 * alpha * pow(l, 2)  + sum) -
                             log(1 + sum / (2 * alpha * pow(l, 2)))) +
                 k * (2 * sum / (2 * alpha * pow(l, 3) + l * sum) -
                      4 * alpha * l * sum / pow(2 * alpha * pow(l, 2) + sum, 2));

  result(1, 0) = result(0, 1);

  double u = 2 * alpha * pow(l, 2) + sum;

  result(1, 1) = k * (sum / u * (1 / alpha - 2 * pow(l, 2) / u) +
                      pow(sum / u - log(1 + sum / (2 * alpha * pow(l, 2))), 2));
  return result;
}


// ****************************************************************************
// * Periodic
// ****************************************************************************


// [[Rcpp::export]]
NumericVector periodicKernel(NumericVector a,
                             NumericVector b,
                             NumericVector hyperParams,
                             List additionalParams) {

  checkHPNames(periodicHPs, hyperParams.names());

  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);

  double l = hyperParams[0];
  double p = hyperParams[1];

  double result = exp(-2 *
                      pow(sin( M_PI *
                                pow(sum, 0.5) /
                                p),
                          2
                      ) /
                      pow(l, 2)
                  );
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector periodicKernelGrad(NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams) {

  checkHPNames(periodicHPs, hyperParams.names());

  // Rcout << a(0) << "; " << b(0) << "; \n";
  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  double l = hyperParams[0];
  double p = hyperParams[1];
  // Rcout << p << "; " << sigma << "; " << l << "; \n";

  double c = M_PI * pow(sum, 0.5) / p;
  double d = 4 * pow(l, -2) * c / p;
  double f = sin(c);
  double g = cos(c);
  double h = exp(-2 * pow(l, -2) * pow(f, 2));

  result["l"] = 4 * pow(f, 2) * pow(l, -3) * h;
  result["p"] = d * g * f * h;
  // Rcout << as<double>(result["sigma"]) << "; " << as<double>(result["l"]) << "; " << as<double>(result["p"]) << "; \n";
  return result;
}

// [[Rcpp::export]]
NumericMatrix periodicKernelHess(NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams) {

  checkHPNames(periodicHPs, hyperParams.names());

  // Rcout << a(0) << "; " << b(0) << "; \n";
  double sumSqrt = pow(sumSQuaredDiffsPartial(a, b, additionalParams), 0.5);

  double k = periodicKernel(a, b, hyperParams, additionalParams)[0];

  NumericMatrix result(hyperParams.size());
  result.attr("dimnames") = List::create(hyperParams.names(), hyperParams.names());

  double l = hyperParams[0];
  double p = hyperParams[1];

  double w = M_PI * sumSqrt / p;
  double u = pow(sin(w), 2);

  double w_p = -M_PI * sumSqrt * pow(p, -2);
  double u_p = w_p * sin(2 * w);

  double w_pp = 2 * M_PI * sumSqrt * pow(p, -3);
  double u_pp = w_pp * sin(2 * w) + 2 * pow(w_p, 2) * cos(2 * w);

  // k_ll
  result(0, 0) = 4 * u * pow(l, -4) * k * (4 * u * pow(l, -2) - 3);

  // k_lp
  result(0, 1) = 4 * pow(l, -3) * u_p * k * (1 - 2 * u * pow(l, -2));
  result(1, 0) = result(0, 1);

  // k_pp
  result(1, 1) = 2 * pow(l, -2) * k * (2 * pow(l, -2) * pow(u_p, 2) - u_pp);

  return result;
}



// ****************************************************************************
// * Constant
// ****************************************************************************


// [[Rcpp::export]]
NumericVector constantKernel(NumericVector a,
                             NumericVector b,
                             NumericVector hyperParams,
                             List additionalParams) {

  checkHPNames(constantHPs, hyperParams.names());

  double result = pow(hyperParams[0], 2);
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector constantKernelGrad(NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams) {

  checkHPNames(constantHPs, hyperParams.names());

  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  result[0] = 2 * hyperParams[0];
  return result;
}

// [[Rcpp::export]]
NumericMatrix constantKernelHess(NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams) {

  checkHPNames(constantHPs, hyperParams.names());

  NumericMatrix result(hyperParams.size());
  result.attr("dimnames") = List::create(hyperParams.names(), hyperParams.names());

  result(0, 0) = 2;

  return result;
}

/*
// ****************************************************************************
// * One-D Linear
// ****************************************************************************


// [[Rcpp::export]]
NumericVector oneDLinearKernel(NumericVector a,
                               NumericVector b,
                               NumericVector hyperParams,
                               List additionalParams) {
  NumericVector dimensionIndex = as<NumericVector>(additionalParams[dimensionIndexName]) - 1;
  double dotProduct = kahanSum((as<NumericVector>(a[dimensionIndex]) - hyperParams["intercept"]) *
                               (as<NumericVector>(b[dimensionIndex]) - hyperParams["intercept"]));
  double result = 1 + pow(hyperParams["sigma_1"], 2) * dotProduct;
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector oneDLinearKernelGrad(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams) {
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);
  double dotProduct = kahanSum((a - hyperParams["intercept"]) * (b - hyperParams["intercept"]));
  double dotProductGrad = kahanSum(2 * hyperParams["intercept"] - a - b);

  result["sigma_1"] = 2 * hyperParams["sigma_1"] * dotProduct;
  result["intercept"] = pow(hyperParams["sigma_1"], 2) * dotProductGrad;
  return result;
}



// ****************************************************************************
// * Generalised Linear
// ****************************************************************************


// [[Rcpp::export]]
NumericVector generalisedLinearKernel(NumericVector a,
                                      NumericVector b,
                                      NumericVector hyperParams,
                                      List additionalParams) {
  NumericVector::iterator weightsIterator = hyperParams.begin();
  weightsIterator++;
  NumericVector weights(weightsIterator, hyperParams.end());
  double result = kahanSum(a * b * pow(weights, 2)) + pow(hyperParams["sigma_0"], 2);
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector generalisedLinearKernelGrad(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams) {
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  result["sigma_0"] = 2 * hyperParams["sigma_0"];
  for (int i=1; i<hyperParams.size(); i++) {
    result[i] = 2 * hyperParams[i] * a[i-1] * b[i-1];
  }
  return result;
}
*/
// ****************************************************************************
// * Generalised Polynomial
// ****************************************************************************


// [[Rcpp::export]]
NumericVector generalisedPolynomialKernel(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams) {

  checkHPNames(generalisedPolynomialHPs, hyperParams.names());

  double degree = additionalParams["degree"];
  double c = hyperParams[0];
  double l = hyperParams[1];

  double result = pow(kahanSum(a * b) / pow(l, 2) + pow(c, 2), degree);
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector generalisedPolynomialKernelGrad(NumericVector a,
                                              NumericVector b,
                                              NumericVector hyperParams,
                                              List additionalParams) {

  checkHPNames(generalisedPolynomialHPs, hyperParams.names());

  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  double degree = additionalParams["degree"];
  double c = hyperParams[0];
  double l = hyperParams[1];

  result[0] = 2 * c;
  result[1] = -2 * kahanSum(a * b) / pow(l, 3);
  if (degree!=1) {
    double innerDeriv = degree * pow(kahanSum(a * b) / pow(l, 2) + pow(c, 2), degree - 1);
    result[0] = result[0] * innerDeriv;
    result[1] = result[1] * innerDeriv;
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix generalisedPolynomialKernelHess(NumericVector a,
                                              NumericVector b,
                                              NumericVector hyperParams,
                                              List additionalParams) {

  checkHPNames(generalisedPolynomialHPs, hyperParams.names());

  NumericMatrix result(hyperParams.size());
  result.attr("dimnames") = List::create(hyperParams.names(), hyperParams.names());

  double d = additionalParams["degree"];
  double c = hyperParams[0];
  double l = hyperParams[1];

  double dot = kahanSum(a * b);
  double inner = dot * pow(l, -2) + pow(c, 2);


  result(1, 1) = 6 * dot * pow(l, -4) * d * pow(inner, d - 1) +
                 4 * pow(dot, 2) * pow(l, -6) * d * (d - 1) * pow(inner, d - 2);
  result(1, 0) = -4 * c * d * dot * pow(l, -3) * (d - 1) * pow(inner, d - 2);
  result(0, 1) = result(0, 1);
  result(0, 0) = 2 * d * pow(inner, d - 1) + 4 * pow(c, 2) * d * (d - 1) * pow(inner, d - 2);

  return result;
}

// ****************************************************************************
// * Polynomial
// ****************************************************************************


// [[Rcpp::export]]
NumericVector polynomialKernel(NumericVector a,
                               NumericVector b,
                               NumericVector hyperParams,
                               List additionalParams) {

  checkHPNames(polynomialHPs, hyperParams.names());

  NumericVector polyHyperParams = NumericVector::create(_["c"] = hyperParams[0],
                                                        _["l"] = 1);

  return generalisedPolynomialKernel(a, b, polyHyperParams, additionalParams);
}

// [[Rcpp::export]]
NumericVector polynomialKernelGrad(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams) {

  checkHPNames(polynomialHPs, hyperParams.names());

  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  NumericVector polyHyperParams = NumericVector::create(_["c"] = hyperParams[0],
                                                        _["l"] = 1);

  NumericVector genGrad = generalisedPolynomialKernelGrad(a, b, polyHyperParams, additionalParams);

  result[0] = genGrad[0];
  return result;
}

// [[Rcpp::export]]
NumericMatrix polynomialKernelHess(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams) {

  checkHPNames(polynomialHPs, hyperParams.names());

  NumericMatrix result(hyperParams.size());

  NumericVector polyHyperParams = NumericVector::create(_["c"] = hyperParams[0],
                                                        _["l"] = 1);

  NumericMatrix genHess = generalisedPolynomialKernelHess(a, b, polyHyperParams, additionalParams);

  result(0, 0) = genHess(0, 0);
  return result;
}

// ****************************************************************************
// * Homogeneous Polynomial
// ****************************************************************************


// [[Rcpp::export]]
NumericVector homogeneousPolynomialKernel(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams) {

  CharacterVector hpNames;
  if (hyperParams.size() == 0) {
    hpNames = CharacterVector::create();
  } else {
    hpNames = hyperParams.names();
  }
  checkHPNames(homogeneousPolynomialHPs, hpNames);

  NumericVector polyHyperParams = NumericVector::create(_["c"] = 0);

  return polynomialKernel(a, b, polyHyperParams, additionalParams);
}

// [[Rcpp::export]]
NumericVector homogeneousPolynomialKernelGrad(NumericVector a,
                                              NumericVector b,
                                              NumericVector hyperParams,
                                              List additionalParams) {

  CharacterVector hpNames;
  if (hyperParams.size() == 0) {
    hpNames = CharacterVector::create();
  } else {
    hpNames = hyperParams.names();
  }
  checkHPNames(homogeneousPolynomialHPs, hpNames);

  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  return result;
}

// [[Rcpp::export]]
NumericMatrix homogeneousPolynomialKernelHess(NumericVector a,
                                              NumericVector b,
                                              NumericVector hyperParams,
                                              List additionalParams) {

  CharacterVector hpNames;
  if (hyperParams.size() == 0) {
    hpNames = CharacterVector::create();
  } else {
    hpNames = hyperParams.names();
  }
  checkHPNames(homogeneousPolynomialHPs, hpNames);

  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericMatrix result(hyperParams.size());

  return result;
}
/*
// ****************************************************************************
// * Changepoint
// ****************************************************************************

double sigmoid(double x,
               double changepoint,
               double transitionRate) {
  return(1 / (1 + exp(transitionRate * (x - changepoint))));
}

// [[Rcpp::export]]
NumericVector changepointKernel(NumericVector a,
                                       NumericVector b,
                                       NumericVector hyperParams,
                                       List additionalParams) {
  NumericVector dimensionIndex = as<NumericVector>(additionalParams[dimensionIndexName]) - 1;
  double changepoint = hyperParams["changepoint"];
  double transitionRate = hyperParams["transitionRate"];

  double result = sigmoid(as<double>(a[dimensionIndex]), changepoint, transitionRate) * sigmoid(as<double>(b[dimensionIndex]), changepoint, transitionRate);
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector changepointKernelGrad(NumericVector a,
                                           NumericVector b,
                                           NumericVector hyperParams,
                                           List additionalParams) {
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);
  NumericVector dimensionIndex = as<NumericVector>(additionalParams[dimensionIndexName]) - 1;
  double c = hyperParams["changepoint"];
  double t = hyperParams["transitionRate"];

  double aval = as<double>(a[dimensionIndex]);
  double bval = as<double>(b[dimensionIndex]);

  double sigmoidvala = sigmoid(aval, c, t);
  double sigmoidvalb = sigmoid(bval, c, t);

  result["changepoint"] = -t * exp(t * (aval - c)) * pow(sigmoidvala, 2) * sigmoidvalb -
                              t * exp(t * (bval - c)) * pow(sigmoidvalb, 2) * sigmoidvala;
  result["transitionRate"] = pow(sigmoidvala, 2) * (aval - c) * exp(t * (aval - c)) * sigmoidvalb +
                              pow(sigmoidvalb, 2) * (bval - c) * exp(t * (bval - c)) * sigmoidvala;
  return result;
}

*/
// ****************************************************************************
// * Random Forest
// ****************************************************************************

// [[Rcpp::export]]
NumericVector randomForestKernel(NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams) {
  double result;
  double matchCount = 0;
  double aCategory;
  double bCategory;

  const char forestName[] = "forest";
  const char heightsName[] = "heights";
  const char isFactorName[] = "isFactor";

  List forest = additionalParams[forestName];
  NumericVector heights = additionalParams[heightsName];
  IntegerVector isFactor = additionalParams[isFactorName];
  double numTrees = forest.length();

  for (int i=0; i<numTrees; i++) {
    NumericVector height = NumericVector::create(heights[i]);
    //Rcout << "Height: " << heights[i] << "\n";
    aCategory = navigateRFTree(forest[i], a, isFactor, height)[0];
    bCategory = navigateRFTree(forest[i], b, isFactor, height)[0];
    //Rcout << "A: " << aCategory << " B: " << bCategory << "\n";
    if (aCategory == bCategory) {
      matchCount++;
    }
  }
  //Rcout << matchCount << "\n";
  result = matchCount / numTrees;
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector randomForestKernelGrad(NumericVector a,
                                     NumericVector b,
                                     NumericVector hyperParams,
                                     List additionalParams) {
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);
  return result;
}

// [[Rcpp::export]]
NumericMatrix randomForestKernelHess(NumericVector a,
                                     NumericVector b,
                                     NumericVector hyperParams,
                                     List additionalParams) {

  CharacterVector hpNames;
  if (hyperParams.size() == 0) {
    hpNames = CharacterVector::create();
  } else {
    hpNames = hyperParams.names();
  }
  checkHPNames(randomForestHPs, hpNames);

  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericMatrix result(hyperParams.size());

  return result;
}

// ****************************************************************************
// * Neural Network
// ****************************************************************************

// [[Rcpp::export]]
NumericVector neuralNetworkKernel(NumericVector a,
                                  NumericVector b,
                                  NumericVector hyperParams,
                                  List additionalParams) {

  double sigma0 = hyperParams["sigma_0"];
  double sigma = hyperParams["sigma"];

  double numerator = 2 * (pow(sigma0, 2) + sum(a * pow(sigma, 2) * b));
  double denominator = sqrt((1 + 2 * (pow(sigma0, 2) + sum(a * pow(sigma, 2) * a))) *
                            (1 + 2 * (pow(sigma0, 2) + sum(b * pow(sigma, 2) * b))));
  double result = 2 / M_PI * asin(numerator / denominator);

  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector neuralNetworkKernelGrad(NumericVector a,
                                      NumericVector b,
                                      NumericVector hyperParams,
                                      List additionalParams) {
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);
  double sigma0 = hyperParams["sigma_0"];
  double sigma = hyperParams["sigma"];
  //Rcout << "sigma0 " << sigma0 << " sigma " << sigma << "\n";
  double aSum = kahanSum(pow(a, 2));
  double bSum = kahanSum(pow(b, 2));
  double abSum = kahanSum(a * b);

  //Rcout << "aSum " << aSum << " bSum " << bSum << " abSum " << abSum << "\n";

  double sigma0_2 = pow(sigma0, 2);
  double sigma_2 = pow(sigma, 2);

  //Rcout << "sigma0_2 " << sigma0_2 << " sigma_2 " << sigma_2 << "\n";

  double subdenominator1 = (1 + 2 * (sigma0_2 + sigma_2 * aSum)) *
    (1 + 2 * (sigma0_2 + sigma_2 * bSum));

  double numerator1 = 4 * sigma0 / sqrt(subdenominator1);

  double subnumerator1 = (sigma0_2 + sigma_2 * abSum) *
                         ( 4 * sigma0 * (1 + 2 * (sigma0_2 + sigma_2 * aSum)) +
                           4 * sigma0 * (1 + 2 * (sigma0_2 + sigma_2 * bSum)));

  double numerator2 = subnumerator1 / sqrt(pow(subdenominator1, 3));
  //Rcout << "subnumerator1 " << subnumerator1 << " pow(subdenominator1, 3/2) " << pow(subdenominator1, 3/2) << "\n";

  double numerator_s0 = 2 * (numerator1 - numerator2);

  double denominator = M_PI * sqrt(1 - 4 * pow(sigma0_2 + sigma_2 * abSum, 2) /
                                     subdenominator1);

  result["sigma_0"] = numerator_s0 / denominator;


  double numerator3 = 4 * sigma * abSum / sqrt(subdenominator1);

  double subnumerator2 = (sigma0_2 + sigma_2 * abSum) *
                          ( 4 * sigma * bSum * (1 + 2 * (sigma0_2 + sigma_2 * aSum)) +
                            4 * sigma * aSum * (1 + 2 * (sigma0_2 + sigma_2 * bSum)));

  double numerator4 = subnumerator2 / sqrt(pow(subdenominator1, 3));

  double numerator_s = 2 * (numerator3 - numerator4);

  result["sigma"] = numerator_s / denominator;

  return result;
}

// ****************************************************************************
// * General Neural Network
// ****************************************************************************

// [[Rcpp::export]]
NumericVector generalNeuralNetworkKernel(NumericVector a,
                                         NumericVector b,
                                         NumericVector hyperParams,
                                         List additionalParams) {

  double dimensions = a.length();
  NumericVector aTilde(dimensions + 1, 1);
  NumericVector bTilde(dimensions + 1, 1);
  for (int i=0; i<dimensions; i++) {
    aTilde[i] = a[i];
    bTilde[i] = b[i];
  }

  double numerator = 2 * (sum(aTilde * pow(hyperParams, 2) * bTilde));
  double denominator = sqrt((1 + 2 * sum(pow(aTilde, 2) * pow(hyperParams, 2))) *
                            (1 + 2 * sum(pow(bTilde, 2) * pow(hyperParams, 2))));
  double result = 2 / M_PI * asin(numerator / denominator);

  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector generalNeuralNetworkKernelGrad(NumericVector a,
                                             NumericVector b,
                                             NumericVector hyperParams,
                                             List additionalParams) {
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  double dimensions = a.length();
  NumericVector aTilde(dimensions + 1, 1);
  NumericVector bTilde(dimensions + 1, 1);
  for (int i=0; i<dimensions; i++) {
    aTilde[i] = a[i];
    bTilde[i] = b[i];
  }

  double num = 2 * (sum(aTilde * pow(hyperParams, 2) * bTilde));

  double denom1 = 1 + 2 * sum(pow(aTilde, 2) * pow(hyperParams, 2));
  double denom2 = 1 + 2 * sum(pow(bTilde, 2) * pow(hyperParams, 2));

  double denom = sqrt(denom1 * denom2);

  double u = num / denom;

  double v = pow(1 - pow(u, 2), -0.5);

  double q = pow(denom, -1);

  for (int i=0; i<=dimensions; i++) {
    double du_dsigma_i = 2 * hyperParams[i] * q * (
      2 * aTilde[i] * bTilde[i] -
        num * ( pow(aTilde[i], 2) / denom1 +
                pow(bTilde[i], 2) / denom2
        )
      );
    result[i] = 2 / M_PI * v * du_dsigma_i;
  }

  return result;
}

// ****************************************************************************
// * Kernel Selector
// ****************************************************************************

// This function selects between kernels. Remember to add new kernels in here.
kernPtr selectKernel(std::string kernelName, bool returnGrad) {

  if (kernelName == "squaredExponential") {
    if (returnGrad) {
      return(squaredExponentialKernelGrad);
    } else {
      return(squaredExponentialKernel);
    }
  } else if (kernelName == "ARD") {
    if (returnGrad) {
      return(ARDKernelGrad);
    } else {
      return(ARDKernel);
    }
  } else if (kernelName == "inverseARD") {
    if (returnGrad) {
      return(inverseARDKernelGrad);
    } else {
      return(inverseARDKernel);
    }
  } else if (kernelName == "rationalQuadratic") {
    if (returnGrad) {
      return(rationalQuadraticKernelGrad);
    } else {
      return(rationalQuadraticKernel);
    }
  } else if (kernelName == "periodic") {
    if (returnGrad) {
      return(periodicKernelGrad);
    } else {
      return(periodicKernel);
    }
  } else if (kernelName == "constant") {
    if (returnGrad) {
      return(constantKernelGrad);
    } else {
      return(constantKernel);
    } /*
  } else if (kernelName == "generalisedLinear") {
    if (returnGrad) {
      return(generalisedLinearKernelGrad);
    } else {
      return(generalisedLinearKernel);
    }
  } else if (kernelName == "oneDLinear") {
    if (returnGrad) {
      return(oneDLinearKernelGrad);
    } else {
      return(oneDLinearKernel);
    }
  } else if (kernelName == "changepoint") {
    if (returnGrad) {
      return(changepointKernelGrad);
    } else {
      return(changepointKernel);
    } */
  } else if (kernelName == "randomForest") {
    if (returnGrad) {
      return(randomForestKernelGrad);
    } else {
      return(randomForestKernel);
    }
  } else if (kernelName == "neuralNetwork") {
    if (returnGrad) {
      return(neuralNetworkKernelGrad);
    } else {
      return(neuralNetworkKernel);
    }
  } else if (kernelName == "generalisedNeuralNetwork") {
    if (returnGrad) {
      return(generalNeuralNetworkKernelGrad);
    } else {
      return(generalNeuralNetworkKernel);
    }
  } else if (kernelName == "generalisedPolynomial") {
    if (returnGrad) {
      return(generalisedPolynomialKernelGrad);
    } else {
      return(generalisedPolynomialKernel);
    }
  } else if (kernelName == "polynomial") {
    if (returnGrad) {
      return(polynomialKernelGrad);
    } else {
      return(polynomialKernel);
    }
  } else if (kernelName == "homogeneousPolynomial") {
    if (returnGrad) {
      return(homogeneousPolynomialKernelGrad);
    } else {
      return(homogeneousPolynomialKernel);
    }
  } else {
    throw std::range_error("Incorrect kernel specified");
  }
}

// This function selects between kernel hessians. Remember to add new kernels in here.
kernHessPtr selectKernelHess(std::string kernelName) {
  kernHessPtr dummyFun;
  if (kernelName == "squaredExponential") {
    return(squaredExponentialKernelHess);
  } else if (kernelName == "ARD") {
    return(ARDKernelHess);
  } else if (kernelName == "inverseARD") {
    return(inverseARDKernelHess);
  } else if (kernelName == "rationalQuadratic") {
    return(rationalQuadraticKernelHess);
  } else if (kernelName == "periodic") {
    return(periodicKernelHess);
  } else if (kernelName == "constant") {
    return(constantKernelHess); /*
  } else if (kernelName == "generalisedLinear") {
    return(dummyFun);
  } else if (kernelName == "oneDLinear") {
    return(dummyFun);
  } else if (kernelName == "changepoint") {
    return(dummyFun); */
  } else if (kernelName == "randomForest") {
    return(randomForestKernelHess);
  } else if (kernelName == "neuralNetwork") {
    return(dummyFun);
  } else if (kernelName == "generalisedNeuralNetwork") {
    return(dummyFun);
  } else if (kernelName == "generalisedPolynomial") {
    return(generalisedPolynomialKernelHess);
  } else if (kernelName == "polynomial") {
    return(polynomialKernelHess);
  } else if (kernelName == "homogeneousPolynomial") {
    return(homogeneousPolynomialKernelHess);
  } else {
    throw std::range_error("Incorrect kernel specified");
  }
}

