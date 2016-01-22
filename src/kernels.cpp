#define _USE_MATH_DEFINES
#include "utils.hpp"
#include "kernels.hpp"
#include <Rcpp.h>
#include "rfTreeUtils.hpp"
#include <math.h>
#include <string>
using namespace Rcpp;


// ****************************************************************************
// * Squared Exponential
// ****************************************************************************

// [[Rcpp::export]]
NumericVector squaredExponentialKernel(NumericVector a,
                                       NumericVector b,
                                       NumericVector hyperParams,
                                       List additionalParams) {
  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);
  double result = exp(-sum / ( 2 * pow(hyperParams["l"], 2)));
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector squaredExponentialKernelGrad(NumericVector a,
                                           NumericVector b,
                                           NumericVector hyperParams,
                                           List additionalParams) {
  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);
  result["l"] = sum / pow(hyperParams["l"], 3) *
    exp(-sum/(2 * pow(hyperParams["l"], 2)));
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
  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);
  double l = hyperParams["l"];
  double alpha = hyperParams["alpha"];
  double result = pow(1 + sum / ( 2 * alpha * pow(l, 2)), -alpha);
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector rationalQuadraticKernelGrad(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams) {
  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  double l = hyperParams["l"];
  double alpha = hyperParams["alpha"];
  //Rcout << sigma << "; " << l << "; " << alpha << ";\n";

  double u = sum / (2 * pow(l, 2));
  double v = 1 + u / alpha;
  double w = sum / pow(l, 3);
  double x = u / (alpha + u);
  //Rcout << u << "; " << v << "; " << w << "; " << x << ";\n";

  result["l"] = w * pow(v, -(alpha + 1));
  result["alpha"] = pow(v, -alpha) * (x - log(v));

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
  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);
  double result = exp(-2 *
                      pow(sin( M_PI *
                                pow(sum, 0.5) /
                                hyperParams["p"]),
                          2
                      ) /
                      pow(hyperParams["l"], 2)
                  );
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector periodicKernelGrad(NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams) {
  // Rcout << a(0) << "; " << b(0) << "; \n";
  double sum = sumSQuaredDiffsPartial(a, b, additionalParams);
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  double p = hyperParams["p"];
  double l = hyperParams["l"];
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



// ****************************************************************************
// * Constant
// ****************************************************************************


// [[Rcpp::export]]
NumericVector constantKernel(NumericVector a,
                             NumericVector b,
                             NumericVector hyperParams,
                             List additionalParams) {
  double result = pow(hyperParams["sigma_0"], 2);
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector constantKernelGrad(NumericVector a,
                                 NumericVector b,
                                 NumericVector hyperParams,
                                 List additionalParams) {
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  result["sigma_0"] = 2 * hyperParams["sigma_0"];
  return result;
}


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

// ****************************************************************************
// * Generalised Polynomial
// ****************************************************************************


// [[Rcpp::export]]
NumericVector generalisedPolynomialKernel(NumericVector a,
                                          NumericVector b,
                                          NumericVector hyperParams,
                                          List additionalParams) {
  double degree = additionalParams["degree"];
  double l = hyperParams["l"];
  double c = hyperParams["c"];

  double result = pow(kahanSum(a * b) / pow(l, 2) + pow(c, 2), degree);
  return NumericVector::create(result);
}

// [[Rcpp::export]]
NumericVector generalisedPolynomialKernelGrad(NumericVector a,
                                              NumericVector b,
                                              NumericVector hyperParams,
                                              List additionalParams) {
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  double degree = additionalParams["degree"];
  double l = hyperParams["l"];
  double c = hyperParams["c"];

  result["c"] = 2 * c;
  result["l"] = -2 * kahanSum(a * b) / pow(l, 3);
  if (degree!=1) {
    double innerDeriv = degree * pow(kahanSum(a * b) / pow(l, 2) + pow(c, 2), degree - 1);
    result["c"] = result["c"] * innerDeriv;
    result["l"] = result["l"] * innerDeriv;
  }
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
  NumericVector polyHyperParams = clone<NumericVector>(hyperParams);
  polyHyperParams["l"] = 1;

  return generalisedPolynomialKernel(a, b, polyHyperParams, additionalParams);
}

// [[Rcpp::export]]
NumericVector polynomialKernelGrad(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams) {
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  NumericVector polyHyperParams = clone<NumericVector>(hyperParams);
  polyHyperParams["l"] = 1;

  NumericVector genGrad = generalisedPolynomialKernelGrad(a, b, polyHyperParams, additionalParams);

  result["c"] = genGrad["c"];
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
  NumericVector polyHyperParams = clone<NumericVector>(hyperParams);
  polyHyperParams["c"] = 0;

  return polynomialKernel(a, b, polyHyperParams, additionalParams);
}

// [[Rcpp::export]]
NumericVector homogeneousPolynomialKernelGrad(NumericVector a,
                                   NumericVector b,
                                   NumericVector hyperParams,
                                   List additionalParams) {
  // Copy hyperParams for our output vector to ensure we return things in the right order
  NumericVector result = clone<NumericVector>(hyperParams);

  return result;
}

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
  const char cacheName[] = "cache";

  bool caching = false;

  List forest = additionalParams[forestName];
  NumericVector heights = additionalParams[heightsName];
  IntegerVector isFactor = additionalParams[isFactorName];
  double numTrees = forest.length();

  if (additionalParams.containsElementNamed(cacheName)) {
    caching = true;
    List cache = additionalParams[cacheName];
  }

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
    }
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
    }
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
