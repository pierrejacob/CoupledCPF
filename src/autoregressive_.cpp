#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;


// [[Rcpp::export]]
NumericMatrix ar_rinit_rcpp(int nparticles, NumericVector & theta, NumericVector & rand, int dimension){
  NumericMatrix xparticles(nparticles, dimension);
  for (int istate = 0; istate < dimension; istate ++){
    for (int iparticle = 0; iparticle < nparticles; iparticle ++){
      xparticles(iparticle, istate) =  rand((istate * nparticles) + iparticle);
    }
  }
  return xparticles;
}

// [[Rcpp::export]]
NumericMatrix ar_rtransition_rcpp(NumericMatrix & xparticles, NumericVector & theta,
                                     int time, NumericVector & rand, int dimension, NumericMatrix & A){
  // it's a bit wasteful to do the construction of A at every transition
  // Eigen::MatrixXd A(dimension, dimension);
  // for (int i = 0; i < dimension; i++){
  //   for (int j = 0; j < dimension; j++){
  //     A(i,j) = std::pow(theta(0), 1 + std::abs(i - j));
  //   }
  // }
  // Eigen::Map<Eigen::MatrixXd> xparticle_e(as<Eigen::Map<Eigen::MatrixXd> >(xparticles));
  // Eigen::MatrixXd new_x = xparticle_e * A;
  int nparticles = xparticles.nrow();
  NumericMatrix new_x(nparticles, dimension);
  std::fill(new_x.begin(), new_x.end(), 0.);
  for (int istate = 0; istate < dimension; istate ++){
    for (int iparticle = 0; iparticle < nparticles; iparticle ++){
      for (int j = 0; j < dimension; j ++){
        new_x(iparticle, istate) +=  xparticles(iparticle, j) *  A(j, istate);
      }
      new_x(iparticle, istate) = rand((istate * nparticles) + iparticle);
    }
  }
  return new_x;
}
