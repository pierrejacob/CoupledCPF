#include <RcppEigen.h>
#include "systematic.h"
#include "compute_cost.h"
#include "median.h"
#include "wasserstein_auto.h"


using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerVector transport_given_cpp(NumericMatrix x1, NumericMatrix x2, NumericVector w1, NumericVector w2, NumericVector uniforms,
                            double epsilon, double desired_alpha, IntegerVector ancestors_ref){
  int nparticles = w1.size();
  double unif_systematic1 = uniforms(2*nparticles);
  
  NumericMatrix M = cost_matrix_(x1, x2);
  List wass_results = wasserstein_auto_(w1, w2, M, epsilon * median_rcpp(M), desired_alpha);
  NumericMatrix couplingmatrix = wass_results["transportmatrix"];
  // fix the marginal
  // compute the row and column sums
  NumericVector u_row(nparticles);
  NumericVector u_col(nparticles);
  for (int i = 0; i < nparticles; i++) {
    double total1 = 0;
    double total2 = 0;
    for (int j = 0; j < nparticles; j++) {
      total1 += couplingmatrix(i, j);
      total2 += couplingmatrix(j, i);
    }
    u_row(i) = total1;
    u_col(i) = total2;
  }
  //
  NumericVector nu(nparticles);
  double alpha = 1;
  for (int i = 0; i < nparticles; i ++){
    nu(i) = std::min(w1(i)/u_row(i), w2(i)/u_col(i));
    alpha = std::min(alpha, nu(i));
  }
  // residuals
  NumericVector R2 = (w2 - alpha * u_col) / (1 - alpha);
  IntegerVector ancestors2(nparticles);
  //
  LogicalVector coupled(nparticles);
  int ncoupled = 0;
  for (int i = 0; i < nparticles; i ++){
    if (uniforms(i) < alpha){
      coupled(i) = true;
      ncoupled ++;
    } else {
      coupled(i) = false;
    }
  }
  // IntegerVector ancestors_as_vec, a_R1, a_R2;
  IntegerVector a_R2, acoupled;
  if (ncoupled > 0){
    acoupled = IntegerVector(ncoupled);
    int coupledcounter = 0;
    for (int i = 0; i < nparticles; i ++){
      if (coupled(i)){
        NumericVector coupling_row = couplingmatrix(ancestors_ref(i),_);
        coupling_row = coupling_row / sum(coupling_row);
        IntegerVector a = systematic_resampling_n_(coupling_row, 1, uniforms(nparticles + i));
        acoupled(coupledcounter) = a(0);
        coupledcounter ++;
      }
    }
  }
  if (ncoupled < nparticles){
    // a_R1 = systematic_resampling_n_(R1, nparticles - ncoupled, unif_systematic2);
    a_R2 = systematic_resampling_n_(R2, nparticles - ncoupled, unif_systematic1);
  }
  int coupledcounter = 0;
  int uncoupledcounter = 0;
  for (int i = 0; i < nparticles; i ++){
    if (coupled(i)){
      ancestors2(i) = acoupled(coupledcounter);
      coupledcounter ++;
    } else {
      ancestors2(i) = a_R2(uncoupledcounter);
      uncoupledcounter ++;
    }
  }
  return ancestors2;
}

