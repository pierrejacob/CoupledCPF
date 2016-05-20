#include <RcppEigen.h>
#include "mvnorm.h"
using namespace Rcpp;
using namespace std;
using namespace Eigen;


// [[Rcpp::export]]
NumericVector stochvol_dmeas_(NumericMatrix xparticles, List & theta, NumericVector & observation, int dimension){
  Eigen::MatrixXd cholCy_(as<Eigen::MatrixXd>(theta["cholCy"]));
  int nparticles = xparticles.ncol();
  NumericVector res(nparticles);
  NumericMatrix expparticles(dimension, nparticles);
  fill(res.begin(), res.end(), 0.);
  for (int i = 0; i < nparticles; i++){
    for (int j = 0; j < dimension; j++){
      res(i) = res(i) + xparticles(j,i);
      expparticles(j,i) = exp(-0.5*xparticles(j,i)) * observation(j);      
    }
  }
  NumericVector mean(dimension);
  fill(mean.begin(), mean.end(), 0.);
  res = -0.5*res + dmvnorm_transpose_cholesky(expparticles, mean, cholCy_);;
  return res;
}

// dmeasurement <- function(xparticles, theta, observation, ...) {
//   return(fast_dmvnorm_transpose_cholesky(exp(-0.5*xparticles) * observation, rep(0, dimension), theta$cholCy) - 0.5*apply(xparticles, 2, sum))
// }
