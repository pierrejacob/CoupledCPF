#ifndef KALMAN_H_
#define KALMAN_H_

class LinearGaussian;

#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;
using namespace std;

class Kalman {
  public:
    Kalman();
  virtual ~Kalman();
  void setLinearGaussian(LinearGaussian*);
  void filtering();
  void smoothing();
  // pointer to model
  LinearGaussian* model;
  // results
  ArrayXXd likelihoodresults;
  ArrayXXd filteringresults;
  ArrayXXd predictionresults;
  ArrayXXd smoothingresults;
  unsigned int nobs;
  // quantities to store
  vector<ArrayXd> xFilterMeans;
  vector<ArrayXXd> xFilterVariances;
  vector<ArrayXd> xSmoothMeans;
  vector<ArrayXXd> xSmoothVariances;
  vector<VectorXd> histm;
  vector<MatrixXd> histC;
  vector<VectorXd> hista;
  vector<MatrixXd> histR;
  vector<VectorXd> histf;
  vector<MatrixXd> histQ;
  vector<ArrayXd> hists;
  vector<ArrayXXd> histS;
  // useful functions:
  ArrayXd incrLL;
  ArrayXd getIncrementalLL(void);
  double getLL(void);
  // convenient
  NumericVector get_filtering_mean(int time_step);
  NumericMatrix get_filtering_variance(int time_step);
  NumericVector get_smoothing_mean(int time_step);
  NumericMatrix get_smoothing_variance(int time_step);
};

#endif /* KALMAN_H_ */
