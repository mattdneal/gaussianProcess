#ifndef __RFTREEUTILS__   // if x.h hasn't been included yet...
#define __RFTREEUTILS__   //   #define this so the compiler knows it has been included

#include <Rcpp.h>
using namespace Rcpp;

NumericVector getTreeHeight(NumericMatrix tree);

NumericVector navigateRFTree(NumericMatrix tree, NumericVector data, IntegerVector isFactor, NumericVector inputHeight);

#endif