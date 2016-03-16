#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <RcppEigen.h>
#include "kalman.h"
#include "lineargaussian.h"

using namespace std;
using namespace Eigen;

// [[Rcpp::export]]
double kalman_loglikelihood_(const List & parameters, const NumericMatrix & observations){
  LinearGaussian LGModel;
  LGModel.set_parameters(parameters);
  LGModel.set_observations(observations);
  Kalman k;
  k.setLinearGaussian(&LGModel);
  k.filtering();
  return k.getLL();
}

// [[Rcpp::export]]
NumericMatrix kalman_filtering_means_(const List & parameters, const NumericMatrix & observations){
  LinearGaussian LGModel;
  LGModel.set_parameters(parameters);
  LGModel.set_observations(observations);
  Kalman k;
  k.setLinearGaussian(&LGModel);
  k.filtering();
  NumericVector one_mean = k.get_filtering_mean(0);
  NumericMatrix kf_means(observations.nrow()+1, one_mean.size());
  for (int t = 0; t <= observations.nrow(); t++){
    kf_means(t,_) = k.get_filtering_mean(t);
  }
  return kf_means;
}

// [[Rcpp::export]]
NumericMatrix kalman_smoothing_means_(const List & parameters, const NumericMatrix & observations){
  LinearGaussian LGModel;
  LGModel.set_parameters(parameters);
  LGModel.set_observations(observations);
  Kalman k;
  k.setLinearGaussian(&LGModel);
  k.filtering();
  k.smoothing();
  NumericVector one_mean = k.get_smoothing_mean(0);
  NumericMatrix ks_means(observations.nrow()+1, one_mean.size());
  for (int t = 0; t <= observations.nrow(); t++){
    ks_means(t,_) = k.get_smoothing_mean(t);
  }
  return ks_means;
}







