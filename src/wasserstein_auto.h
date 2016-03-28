#ifndef WASS_H_
#define WASS_H_
#include <RcppEigen.h>
using namespace Rcpp;

List wasserstein_auto_(NumericVector p_, NumericVector q_, NumericMatrix cost_matrix_,
                       double epsilon, double desired_alpha);

#endif


