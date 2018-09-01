#include <RcppEigen.h>
#include "multinomial.h"
using namespace Rcpp;
using namespace std;

inline int randWrapper(const int n) {
  RNGScope scope;
  NumericVector ru = runif(1);
  return floor(ru(0)*n);
}

// [[Rcpp::export]]
IntegerVector multinomial_resampling_n_(const NumericVector & weights, int ndraws){
  RNGScope scope;
  int nparticles = weights.size();
  IntegerVector ancestors(ndraws);
  NumericVector cumsumw = cumsum(weights);
  NumericVector uniforms = runif(ndraws);
  double sumw = cumsumw(nparticles - 1);
  double lnMax = 0;
  int j = nparticles;
  for (int i = ndraws; i > 0; i--){
    lnMax += log(uniforms(i-1)) / i;
    uniforms(i-1) = sumw * exp(lnMax);
    while (j > 0 && uniforms(i-1) < cumsumw(j-1)){
      j --;
    }
    ancestors(i-1) = j;
  }
  std::random_shuffle(ancestors.begin(), ancestors.end(), randWrapper);
  return ancestors;
}

// [[Rcpp::export]]
IntegerMatrix coupled_multinomial_resampling_n_(const NumericVector & weights1, const NumericVector & weights2, int ndraws){
  RNGScope scope;
  int nparticles = weights1.size();
  IntegerVector ancestors1(ndraws);
  IntegerVector ancestors2(ndraws);
  IntegerMatrix ancestors(ndraws, 2);
  NumericVector cumsumw1 = cumsum(weights1);
  NumericVector cumsumw2 = cumsum(weights2);
  NumericVector uniforms = runif(ndraws);
  double lnMax = 0;
  int j1 = nparticles;
  int j2 = nparticles;
  for (int i = ndraws; i > 0; i--){
    lnMax += log(uniforms(i-1)) / i;
    uniforms(i-1) =  exp(lnMax);
    while (j1 > 0 && uniforms(i-1) < cumsumw1(j1-1)){
      j1 --;
    }
    while (j2 > 0 && uniforms(i-1) < cumsumw2(j2-1)){
      j2 --;
    }
    ancestors1(i-1) = j1;
    ancestors2(i-1) = j2;
  }
  
  std::vector<int> indexes;
  indexes.reserve(ndraws);
  for (int i = 0; i < ndraws; ++i){
    indexes.push_back(i);
  }
  std::random_shuffle(indexes.begin(), indexes.end(), randWrapper);
  for (int i = 0; i < ndraws; i++){
    ancestors(i,0) = ancestors1[indexes[i]];
    ancestors(i,1) = ancestors2[indexes[i]];
  }
  return ancestors;
}
