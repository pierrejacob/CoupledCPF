#ifndef MULTINOMIALRES_H_
#define MULTINOMIALRES_H_
#include <RcppEigen.h>
using namespace Rcpp;

IntegerVector multinomial_resampling_n_(const NumericVector & weights, int ndraws);
IntegerMatrix coupled_multinomial_resampling_n_(const NumericVector & weights1, const NumericVector & weights2, int ndraws);

#endif

