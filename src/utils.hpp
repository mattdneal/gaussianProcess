#ifndef __UTILS__   // if x.h hasn't been included yet...
#define __UTILS__   //   #define this so the compiler knows it has been included

#include <Rcpp.h>
using namespace Rcpp;

/* ****************************************************************************
 *   Constants - this needs to be char* not std::string for Rcpp, and 
 *   "const char * const" makes both the char and the pointer constant.
 *****************************************************************************/
 
const char * const dimensionIndicesName = "dimensionIndices";
const char * const dimensionIndexName = "dimensionIndex";

// ****************************************************************************

double kahanSum(NumericVector summands);

double sumSQuaredDiffs(NumericVector a, NumericVector b);

double sumSQuaredDiffsPartial(NumericVector a, 
                              NumericVector b, 
                              List additionalParams);

int bitWiseAnd(NumericVector a, NumericVector b);

int bitWiseAnd(int a, int b);

#endif