#ifndef COMPCOST_H_
#define COMPCOST_H_
#include <RcppEigen.h>
using namespace Rcpp;

NumericMatrix cost_matrix_(const NumericMatrix & x, const NumericMatrix & y);

#endif 

