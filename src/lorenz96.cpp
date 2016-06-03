#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <Rcpp.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array< double , 8> state_type;

class lorenz_ {
  
  double F;
  
  
  public:
    lorenz_( double F) : F(F) { }
  // dx_n/dt = x_n-1 (x_n+1 - x_n-2) - x_n + F
  void operator() ( const state_type &x , state_type &dxdt , const double /* t */ )
  {
    dxdt[0] = x[7] * (x[1] - x[6]) - x[0] + F;
    dxdt[1] = x[0] * (x[2] - x[7]) - x[1] + F;
    dxdt[2] = x[1] * (x[3] - x[0]) - x[2] + F;
    dxdt[3] = x[2] * (x[4] - x[1]) - x[3] + F;
    dxdt[4] = x[3] * (x[5] - x[2]) - x[4] + F;
    dxdt[5] = x[4] * (x[6] - x[3]) - x[5] + F;
    dxdt[6] = x[5] * (x[7] - x[4]) - x[6] + F;
    dxdt[7] = x[6] * (x[0] - x[5]) - x[7] + F;
  }
};

struct push_back_state_and_time
{
  std::vector< state_type >& m_states;
  std::vector< double >& m_times;
  
  push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
  : m_states( states ) , m_times( times ) { }
  
  void operator()( const state_type &x , double t )
  {
    m_states.push_back( x );
    m_times.push_back( t );
  }
};

using namespace Rcpp;

// // [[Rcpp::export]]
// NumericVector one_step_lorenz_(double P, double Z,  double t, double alpha, double c, double e, double ml, double mq ) {
  //   state_type x = { P , Z }; // initial conditions
  //   vector<state_type> x_vec;
  //   vector<double> times;
  //   lorenz_ lorenz_instance(alpha, c, e, ml, mq);
  //   size_t steps = integrate( lorenz_instance , x , t , t + 1, 1.0 , push_back_state_and_time(x_vec, times));
  //   NumericVector result = NumericVector::create(0, 0);
  //   result(0) = x_vec[steps][0];
  //   result(1) = x_vec[steps][1];
  //   return result;
  // }

// [[Rcpp::export]]
NumericMatrix one_step_lorenz_vector(NumericMatrix xparticles, double tstart, double tend, double h, NumericVector parameters){
  double F = parameters[0];
  NumericMatrix result(8, xparticles.cols());
  for (int i = 0; i < xparticles.cols(); i++){
    double x1 = xparticles(0, i);
    double x2 = xparticles(1, i);
    double x3 = xparticles(2, i);
    double x4 = xparticles(3, i);
    double x5 = xparticles(4, i);
    double x6 = xparticles(5, i);
    double x7 = xparticles(6, i);
    double x8 = xparticles(7, i);
    
    state_type x = { x1, x2, x3, x4, x5, x6, x7, x8 }; // initial conditions
    vector<state_type> x_vec;
    vector<double> times;
    lorenz_ lorenz_instance(F);
    size_t steps = integrate( lorenz_instance , x , tstart , tend, h, push_back_state_and_time(x_vec, times));
    result(0, i) = x_vec[steps][0];
    result(1, i) = x_vec[steps][1];
    result(2, i) = x_vec[steps][2];
    result(3, i) = x_vec[steps][3];
    result(4, i) = x_vec[steps][4];
    result(5, i) = x_vec[steps][5];
    result(6, i) = x_vec[steps][6];
    result(7, i) = x_vec[steps][7];
  }
  return result;
}

// [[Rcpp::export]]
NumericVector lorenz_generate_randomness_cpp(int nparticles, int datalength){
  RNGScope scope;
  NumericVector normal_draws = rnorm((8 + 8 * datalength) * nparticles, 0, 1);
  return normal_draws;
}
// 
// [[Rcpp::export]]
NumericVector lorenz_perturb_randomness_cpp(const NumericVector & randomness, double rho){
  RNGScope scope;
  int l = randomness.size();
  NumericVector newrand(l);
  double v = sqrt(1.0 - rho*rho);
  NumericVector normal_draws = rnorm(l, 0, 1);
  for (int i = 0; i < l; i ++){
    newrand(i) =  rho * randomness(i) + v * normal_draws(i);
  }
  return newrand;
}
