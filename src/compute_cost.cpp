#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <RcppEigen.h>
#include "compute_cost.h"
using namespace Rcpp;
using namespace std;
using namespace Eigen;

// [[Rcpp::export]]
NumericMatrix compute_cost1_(const NumericVector & x, const NumericVector & y){
  NumericMatrix cost(x.size(), y.size());
  for (int ix = 0; ix < x.size(); ix ++){
    for (int iy = 0; iy < y.size(); iy ++){
      cost(ix,iy) = std::abs(x(ix)-y(iy));
    }  
  }
  return cost;
}


// [[Rcpp::export]]
NumericMatrix compute_cost2_(const NumericVector & x, const NumericVector & y){
  NumericMatrix cost(x.size(), y.size());
  for (int ix = 0; ix < x.size(); ix ++){
    for (int iy = 0; iy < y.size(); iy ++){
      cost(ix,iy) = std::pow(x(ix)-y(iy), 2);
    }  
  }
  return cost;
}

// [[Rcpp::export]]
NumericMatrix cost_matrix_(const NumericMatrix & x, const NumericMatrix & y){
  NumericMatrix cost(x.rows(), y.rows());
  for (int ix = 0; ix < x.rows(); ix ++){
    for (int iy = 0; iy < y.rows(); iy ++){
      cost(ix,iy) = 0;
      // loop over components, stored in columns
      for (int id = 0; id < x.cols(); id ++){
        // cost(ix,iy) += std::abs(x(ix,id)-y(iy,id));
        cost(ix,iy) += std::pow(x(ix,id)-y(iy,id),2);
      }
      cost(ix,iy) = std::sqrt(cost(ix,iy));
    }  
  }
  return cost;
}

// [[Rcpp::export]]
NumericMatrix square_cost_matrix_(const NumericMatrix & x, const NumericMatrix & y){
  NumericMatrix cost(x.rows(), y.rows());
  for (int ix = 0; ix < x.rows(); ix ++){
    for (int iy = 0; iy < y.rows(); iy ++){
      cost(ix,iy) = 0;
      // loop over components, stored in columns
      for (int id = 0; id < x.cols(); id ++){
        // cost(ix,iy) += std::abs(x(ix,id)-y(iy,id));
        cost(ix,iy) += std::pow(x(ix,id)-y(iy,id),2);
      }
    }  
  }
  return cost;
}
