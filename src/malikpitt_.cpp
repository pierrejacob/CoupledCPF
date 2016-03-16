#include <math.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;

// [[Rcpp::export]]
NumericVector malikpitt_(const NumericVector & xparticles, const NumericVector & weights, double u){
  RNGScope scope;
  int nparticles = weights.size();
  IntegerVector regions(nparticles);
  NumericVector slope(nparticles);
  NumericVector new_x(nparticles);
  u = u / (double) nparticles;
  int j = 0;
  double csw = 0;
  for (int i = 0; i < nparticles; i++){
    csw += weights(i);
    while (u <= csw && j < nparticles){
      regions(j) = i;
      slope(j) = (u - (csw - weights(i))) / weights(i);
      j += 1;
      u += 1 / (double) nparticles;
    }
  }
  for (int i = 0; i < nparticles; i++){
    if (regions(i) == 0 || regions(i) == nparticles - 1){
      new_x(i) = xparticles(i);
    } else {
      new_x(i) = (xparticles(regions(i) + 1) - xparticles(regions(i))) * slope(i) + xparticles(regions(i));
    }
  }
  return new_x;
}
