#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <RcppEigen.h>
#include "kalman.h"
#include "lineargaussian.h"

using namespace std;
using namespace Eigen;

Kalman::Kalman() {
}

Kalman::~Kalman() {
}

void Kalman::setLinearGaussian(LinearGaussian* model){
  this->model = model;
  this->nobs = this->model->observations.rows();
}

void Kalman::filtering(){
  /**
   * * KF computes the Kalman filter for the given observations and the given parameter
   * * Notation taken from Petris, Petrone and Campagnoli: Dynamical Linear Models with R.
   * * (page 41)
   * * Y_t = F_t Theta_t + v_t
   * * Theta_t = G_t Theta_t-1 + w_t
   * * v_t ~ N(0, V_t); w_t ~ N(0, W_t)
   * * Theta_0 ~ N(m_0, C_0)
   * *
   * * Kalman recursion: Proposition 2.2, page 53 for filtering
   * * Kalman smoothing: Proposition 2.4, page 61 for smoothing
   * */
  this->xFilterMeans.clear(); this->xFilterVariances.clear();
  histm.clear(); histC.clear(); hista.clear(); histR.clear();
  histf.clear(); histQ.clear();
  VectorXd m(model->XDimension);
  MatrixXd C(model->XDimension, model->XDimension);
  VectorXd a(model->XDimension);
  MatrixXd R(model->XDimension, model->XDimension);
  VectorXd f(model->YDimension);
  MatrixXd Q(model->YDimension, model->YDimension);
  m = model->m_0;
  C = model->C_0;
  this->xFilterMeans.push_back(m.array());
  this->xFilterVariances.push_back(C.array());
  histm.push_back(m); histC.push_back(C);
  VectorXd e(model->YDimension);
  // incremental likelihoods
  incrLL = ArrayXd::Zero(this->nobs);
  for (unsigned int t = 0; t < this->nobs; t ++){
    a = model->G * m;
    R = model->G * C * model->G.transpose() + model->W;
    hista.push_back(a); histR.push_back(R);
    f = model->F * a;
    Q = model->F * R * model->F.transpose() + model->V;
    histf.push_back(f); histQ.push_back(Q);
    e = this->model->observations.row(t);
    e = e - f;
    m = a + R * (model->F).transpose() * Q.inverse() * e;
    C = R - R * (model->F).transpose() * Q.inverse() * (model->F) * R;
    histm.push_back(m); histC.push_back(C);
    this->xFilterMeans.push_back(m.array());
    this->xFilterVariances.push_back(C.array());
    incrLL(t) = - model->YDimension / 2. * log(2 * M_PI) - 0.5 * log(Q.determinant())
      - 0.5 * e.transpose() * Q.inverse() * e;
  }
}

void Kalman::smoothing(){
  // this function assumes that the filtering function has been called already
  hists.clear(); histS.clear();
  // quantities introduced in Proposition 2.4 for smoothing
  VectorXd s(model->XDimension);
  MatrixXd S(model->XDimension, model->XDimension);
  s = this->histm[this->nobs];
  S = this->histC[this->nobs];
  hists.push_back(s.array());
  histS.push_back(S.array());
  MatrixXd tmpinvR;
  for (int t = (this->nobs - 1); t >= 0; t --){
    tmpinvR = histR[t].inverse();
    s = histm[t] + histC[t] * (model->G).transpose() * tmpinvR * (s - hista[t]);
    S = histC[t] - histC[t] * (model->G).transpose() * tmpinvR * (histR[t] - S) * tmpinvR * model->G * histC[t];
    hists.push_back(s.array());
    histS.push_back(S.array());
  }
  vector<ArrayXd> reversehists = this->hists;
  reverse(reversehists.begin(), reversehists.end());
  xSmoothMeans = reversehists;
  vector<ArrayXXd> reversehistS = this->histS;
  reverse(reversehistS.begin(), reversehistS.end());
  xSmoothVariances = reversehistS;
}


ArrayXd Kalman::getIncrementalLL(void){
  return this->incrLL;
}

double Kalman::getLL(void){
  return this->getIncrementalLL().sum();
}

NumericVector Kalman::get_filtering_mean(int time_step){
  return wrap(this->xFilterMeans[time_step]);
}

NumericMatrix Kalman::get_filtering_variance(int time_step){
  return wrap(this->xFilterVariances[time_step]);
}

NumericVector Kalman::get_smoothing_mean(int time_step){
  return wrap(this->xSmoothMeans[time_step]);
}

NumericMatrix Kalman::get_smoothing_variance(int time_step){
  return wrap(this->xSmoothVariances[time_step]);
}

NumericMatrix get_F(LinearGaussian* LinearGaussian){return wrap(LinearGaussian->F);}
NumericMatrix get_G(LinearGaussian* LinearGaussian){return wrap(LinearGaussian->G);}
NumericMatrix get_V(LinearGaussian* LinearGaussian){return wrap(LinearGaussian->V);}
NumericMatrix get_W(LinearGaussian* LinearGaussian){return wrap(LinearGaussian->W);}
int get_nobs(LinearGaussian* LinearGaussian){return LinearGaussian->observations.rows();}
NumericVector get_incremental_ll(Kalman* kalman){return wrap(kalman->getIncrementalLL());}


RCPP_EXPOSED_CLASS(LinearGaussian)
  RCPP_EXPOSED_CLASS(Kalman)

  RCPP_MODULE(kalman_mod) {
    class_<LinearGaussian>( "LinearGaussian" )
    .constructor()
    .method( "set_F", &LinearGaussian::set_F)
    .method( "set_G", &LinearGaussian::set_G)
    .method( "set_V", &LinearGaussian::set_V)
    .method( "set_W", &LinearGaussian::set_W)
    .method( "get_F", &get_F)
    .method( "get_G", &get_G)
    .method( "get_V", &get_V)
    .method( "get_W", &get_W)
    .method( "setLinearGaussianMatrices", &LinearGaussian::setLinearGaussianMatrices)
    .method( "set_parameters", &LinearGaussian::set_parameters)
    .method( "set_multivariate_parameters", &LinearGaussian::set_multivariate_parameters)
    .method( "set_observations", &LinearGaussian::set_observations)
    .method( "get_nobs", &get_nobs)
    ;

    class_<Kalman>( "Kalman" )
      .constructor()
      .method( "setLinearGaussian", &Kalman::setLinearGaussian)
      .method( "filtering", &Kalman::filtering)
      .method( "smoothing", &Kalman::smoothing)
      .method( "get_filtering_mean", &Kalman::get_filtering_mean)
      .method( "get_filtering_variance", &Kalman::get_filtering_variance)
      .method( "get_smoothing_mean", &Kalman::get_smoothing_mean)
      .method( "get_smoothing_variance", &Kalman::get_smoothing_variance)
      .method( "getLL", &Kalman::getLL)
      .method( "get_incremental_ll", &get_incremental_ll)
      .field( "nobs", &Kalman::nobs)
    ;

  }
