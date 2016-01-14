#include <Rcpp.h>
#include "utils.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
double kahanSum(NumericVector summands) {
  double sum = 0.0;
  double c = 0.0;
  double y;
  double t;
  for (NumericVector::iterator it = summands.begin(); it != summands.end(); ++it) {
    //sum = sum + *it;
    y = *it - c;
    t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
  return sum;
}


// [[Rcpp::export]]
double sumSQuaredDiffs(NumericVector a, NumericVector b) {
  NumericVector squaredDiffs = pow(a - b, 2);
  return(kahanSum(squaredDiffs));
}

// [[Rcpp::export]]
double sumSQuaredDiffsPartial(NumericVector a, 
                              NumericVector b, 
                              List additionalParams) {
  double sum;
  if (additionalParams.containsElementNamed(dimensionIndicesName)) {
    NumericVector dimensionIndices = as<NumericVector>(additionalParams[dimensionIndicesName]);
    // The interface is in R so indices start counting from 1. Shift them down by 1 to make them C++-ey
    dimensionIndices = dimensionIndices - 1;
    sum = sumSQuaredDiffs(a[dimensionIndices], b[dimensionIndices]);
  } else {
    sum = sumSQuaredDiffs(a, b);
  }
  return(sum);
}

// [[Rcpp::export]]
int bitWiseAnd(NumericVector a, NumericVector b) {
  NumericVector ind(1);
  
  unsigned int a1 = as<unsigned int>(a[ind]);
  unsigned int b1 = as<unsigned int>(b[ind]);
  int c = 0;
  c = a1 & b1;
  
  return(c);
}

int bitWiseAnd(int a, int b) {
  
  unsigned int a1 = a;
  unsigned int b1 = b;
  int c = 0;
  c = a1 & b1;
  
  return(c);
}
