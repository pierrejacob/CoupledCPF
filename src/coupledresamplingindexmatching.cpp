#include <RcppEigen.h>
#include "multinomial.h"

using namespace Rcpp;
using namespace std;

// Something's off in this commented code, but I'm not sure why;
// the code in 'coupledresampling-indexmatching.R' is basically the same
// but appears to give different, and more satisfactory results

// // [[Rcpp::export]]
// IntegerMatrix indexmatching_cpp(const NumericVector & w1, const NumericVector & w2, int ndraws){
//   RNGScope scope;
//   int nparticles = w1.size();
//   NumericVector uniforms = runif(ndraws);
//   NumericVector nu(nparticles);
//   double alpha = 0.;
//   for (int i = 0; i < nparticles; i ++){
//     if (w1(i) < w2(i)){
//       nu(i) = w1(i);
//     } else {
//       nu(i) = w2(i);
//     }
//     alpha += nu(i);
//   }
//   // minorization proba
//   NumericVector mu = nu / alpha;
//   // residuals
//   NumericVector R1 = (w1 - nu) / (1. - alpha);
//   NumericVector R2 = (w2 - nu) / (1. - alpha);
//   // 
//   IntegerMatrix ancestors(ndraws, 2);
//   //
//   IntegerVector a(1);
//   IntegerMatrix a_R(1,2);
//   for (int i = 0; i < ndraws; i ++){
//     if (uniforms(i) < alpha){
//       a = multinomial_resampling_n_(mu, 1);
//       ancestors(i,0) = a(0);
//       ancestors(i,1) = a(0);
//     } else {
//       a_R = coupled_multinomial_resampling_n_(R1, R2, 1);
//       ancestors(i,0) = a_R(0,0);
//       ancestors(i,1) = a_R(0,1);
//     }
//   }
//   return ancestors;
// }

