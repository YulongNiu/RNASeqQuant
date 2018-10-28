#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec Logistic(const arma::vec& x) {
  return 1 / (1 + exp(-x));
}

// [[Rcpp::export]]
arma::vec Softplus1(const arma::vec& x) {

}
