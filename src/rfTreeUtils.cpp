#include <Rcpp.h>
#include "utils.hpp"
#include <stack>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector getTreeHeight(NumericMatrix tree) {
  const int leftChild = 1;
  const int rightChild = 2;
  const int nullChild = 0;
  
  int currentHeight = 0;
  int maxHeight = 0;
  
  std::stack<int> priorIndexStack;
  std::stack<int> chosenChildStack;
  
  priorIndexStack.push(0);
  chosenChildStack.push(nullChild);
  
  int currentIndex = 1;
  int lastChosenChild = nullChild;
  
  // Stack indicating left or right daughter selected when last visited? Then no need to alter tree.
  while (currentIndex != 0) {
    maxHeight = std::max(maxHeight, currentHeight);
    if (((tree(currentIndex - 1, 0) == 0) & (tree(currentIndex - 1, 1) == 0)) | 
         (lastChosenChild == rightChild) | 
         ((lastChosenChild == leftChild) & (tree(currentIndex - 1, 1) == 0))) {
      // backtracking - pop the parent off the stack
      currentIndex = priorIndexStack.top();
      priorIndexStack.pop();
      
      lastChosenChild = chosenChildStack.top();
      chosenChildStack.pop();
      
      // decrease our height
      currentHeight -= 1;
    } else {
      // Push this level onto the stack
      priorIndexStack.push(currentIndex);
      currentHeight += 1;
      
      // Choose a child node to investigate
      if ((tree(currentIndex - 1, 0) != 0) & (lastChosenChild == nullChild)) {
        currentIndex = tree(currentIndex - 1, 0);
        chosenChildStack.push(leftChild);
      } else  {
        currentIndex = tree(currentIndex - 1, 1);
        chosenChildStack.push(rightChild);
      }
      
      lastChosenChild = nullChild;
    }
  }
  return(NumericVector::create(maxHeight));
}

// [[Rcpp::export]]
NumericVector navigateRFTree(NumericMatrix tree, NumericVector data, IntegerVector isFactor, NumericVector inputHeight) {
  int currentHeight = 0;
  int currentIndex = 1;
  int height = inputHeight[0];
  //Rcout << height << "\n";
  const int leftIndex = 0;
  const int rightIndex = 1;
  const int splitVarIndex = 2;
  const int splitPointIndex = 3;
  
  
  int splitVar = tree(currentIndex - 1, splitVarIndex);
  
  if (height == -1) {
    height = tree.nrow();
  }
  
  while ((currentHeight < height) & (splitVar != 0)) {
    //Rcout << currentHeight << "\n";
    //Rcout << currentIndex << "\n";
    splitVar = tree(currentIndex - 1, splitVarIndex);
    if (splitVar == 0) {
      // We've hit a leaf node - break out
      break;
    }
    
    double dataSplitVar = data[splitVar - 1];
    
    if (isFactor[splitVar - 1]) {
      // Factor
      int factorSplit = tree(currentIndex - 1, splitPointIndex);
      int factorLevel = dataSplitVar;
      
      //Rcout << factorSplit << " "<< factorLevel << "\n";
      //Rcout << bitWiseAnd(factorSplit, factorLevel);
      if (bitWiseAnd(factorSplit, factorLevel) != 0) {
        currentIndex = tree(currentIndex - 1, leftIndex);
        //Rcout << "l\n";
      } else {
        currentIndex = tree(currentIndex - 1, rightIndex);
        //Rcout << "r\n";
      }
      
    } else {
      double splitPoint = tree(currentIndex - 1, splitPointIndex);
      
      if (dataSplitVar <= splitPoint) {
        currentIndex = tree(currentIndex - 1, leftIndex);
      } else {
        currentIndex = tree(currentIndex - 1, rightIndex);
      }
    } 
    
    currentHeight++;
  }
  
  //Rcout << currentIndex << "\n";
  return(NumericVector::create(currentIndex));
}
