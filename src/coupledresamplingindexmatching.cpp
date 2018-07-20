#include <RcppEigen.h>
#include "multinomial.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerMatrix indexmatching_cpp(NumericVector w1, NumericVector w2, int ndraws){
  RNGScope scope;
  int nparticles = w1.size();
  NumericVector uniforms = runif(nparticles);
  NumericVector nu(nparticles);
  double alpha = 0;
  for (int i = 0; i < nparticles; i ++){
    nu(i) = std::min(w1(i), w2(i));
    alpha += nu(i);
  }
  // minorization proba
  NumericVector mu = nu / alpha;
  // residuals
  NumericVector R1 = (w1 - nu) / (1 - alpha);
  NumericVector R2 = (w2 - nu) / (1 - alpha);
  IntegerMatrix ancestors(ndraws, 2);
  //
  LogicalVector coupled(ndraws);
  int ncoupled = 0;
  for (int i = 0; i < ndraws; i ++){
    if (uniforms(i) < alpha){
      coupled(i) = true;
      ncoupled ++;
    } else {
      coupled(i) = false;
    }
  }
  
  IntegerVector a, a_R1, a_R2;
  if (ncoupled > 0){
    a = multinomial_resampling_n_(mu, ncoupled);
  }
  IntegerMatrix a_R;
  if (ncoupled < ndraws){
    a_R = coupled_multinomial_resampling_n_(R1, R2, ndraws - ncoupled);
  }
  int coupledcounter = 0;
  int uncoupledcounter = 0;
  for (int i = 0; i < ndraws; i ++){
    if (coupled(i)){
      ancestors(i, 0) = a(coupledcounter);
      ancestors(i, 1) = a(coupledcounter);
      coupledcounter ++;
    } else {
      ancestors(i, 0) = a_R(uncoupledcounter,0);
      ancestors(i, 1) = a_R(uncoupledcounter,1);
      uncoupledcounter ++;
    }
  }
  return ancestors;
}

// ////// [[Rcpp::export]]
// IntegerMatrix indexmatching_cpp(int ndraws, NumericVector w1, NumericVector w2, NumericVector uniforms){
//   int nparticles = w1.size();
//   double unif_systematic1 = uniforms(ndraws);
//   double unif_systematic2 = uniforms(ndraws + 1);
//   
//   NumericVector nu(nparticles);
//   double alpha = 0;
//   for (int i = 0; i < nparticles; i ++){
//     nu(i) = std::min(w1(i), w2(i));
//     alpha += nu(i);
//   }
//   // minorization proba
//   NumericVector mu = nu / alpha;
//   // residuals
//   NumericVector R1 = (w1 - nu) / (1 - alpha);
//   NumericVector R2 = (w2 - nu) / (1 - alpha);
//   IntegerMatrix ancestors(ndraws, 2);
//   //
//   LogicalVector coupled(ndraws);
//   int ncoupled = 0;
//   for (int i = 0; i < ndraws; i ++){
//     if (uniforms(i) < alpha){
//       coupled(i) = true;
//       ncoupled ++;
//     } else {
//       coupled(i) = false;
//     }
//   }
//   IntegerVector a, a_R1, a_R2;
//   if (ncoupled > 0){
//     a = systematic_resampling_n_(mu, ncoupled, unif_systematic1);
//   }
//   if (ncoupled < ndraws){
//     a_R1 = systematic_resampling_n_(R1, ndraws - ncoupled, unif_systematic2);
//     a_R2 = systematic_resampling_n_(R2, ndraws - ncoupled, unif_systematic2);
//   }
//   int coupledcounter = 0;
//   int uncoupledcounter = 0;
//   for (int i = 0; i < ndraws; i ++){
//     if (coupled(i)){
//       ancestors(i, 0) = a(coupledcounter);
//       ancestors(i, 1) = a(coupledcounter);
//       coupledcounter ++;
//     } else {
//       ancestors(i, 0) = a_R1(uncoupledcounter);
//       ancestors(i, 1) = a_R2(uncoupledcounter);
//       uncoupledcounter ++;
//     }
//   }
//   return ancestors;
// }

