#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix create_A_(double alpha, int d){
  NumericMatrix A(d, d);
  for (int i = 0; i < d; i++){
    for (int j = 0; j < d; j++){
      A(i,j) = std::pow(alpha, 1 + std::abs(i - j));
    }
  }
  return A;
}