#include <RcppEigen.h>
#include "systematic.h"
#include "compute_cost.h"
#include "median.h"
#include "wasserstein_auto.h"


using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerMatrix transport_cpp(NumericMatrix x1, NumericMatrix x2, NumericVector w1, NumericVector w2, NumericVector uniforms,
                            double epsilon, double desired_alpha){
  int nparticles = w1.size();
  double unif_systematic1 = uniforms(nparticles);
  double unif_systematic2 = uniforms(nparticles + 1);
  
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
  // we compute 
  // nu <- pmin(normweights1 / u_1, normweights2 / u_2)
  // alpha <- min(nu)
  // R1 <- (normweights1 - alpha * u_1) / (1 - alpha)
  // R2 <- (normweights2 - alpha * u_2) / (1 - alpha)
  NumericVector nu(nparticles);
  double alpha = 1;
  for (int i = 0; i < nparticles; i ++){
    nu(i) = std::min(w1(i)/u_row(i), w2(i)/u_col(i));
    alpha = std::min(alpha, nu(i));
  }
  // residuals
  NumericVector R1 = (w1 - alpha * u_row) / (1 - alpha);
  NumericVector R2 = (w2 - alpha * u_col) / (1 - alpha);
  IntegerMatrix ancestors(nparticles, 2);
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
  IntegerVector ancestors_as_vec, a_R1, a_R2;
  if (ncoupled > 0){
    // stack coupling matrix in a vector
    NumericVector couplingvec(nparticles*nparticles);
    int irow, icol;
    for (int i = 0; i < nparticles*nparticles; i++){
      irow = i % nparticles; 
      icol = (i - irow) / nparticles;
      couplingvec(i) = couplingmatrix(irow, icol);
    }
    // then sample using standard systematic resampling 
    ancestors_as_vec = systematic_resampling_n_(couplingvec, ncoupled, unif_systematic1);
  }
  if (ncoupled < nparticles){
    a_R1 = systematic_resampling_n_(R1, nparticles - ncoupled, unif_systematic2);
    a_R2 = systematic_resampling_n_(R2, nparticles - ncoupled, unif_systematic2);
  }
  int coupledcounter = 0;
  int uncoupledcounter = 0;
  for (int i = 0; i < nparticles; i ++){
    if (coupled(i)){
      int a = ancestors_as_vec(coupledcounter);
      ancestors(i, 0) = a % nparticles;
      ancestors(i, 1) = (a - ancestors(i, 0)) / nparticles;;
      coupledcounter ++;
    } else {
      ancestors(i, 0) = a_R1(uncoupledcounter);
      ancestors(i, 1) = a_R2(uncoupledcounter);
      uncoupledcounter ++;
    }
  }
  return ancestors;
}

