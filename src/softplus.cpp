#include <RcppArmadillo.h>

#include <cmath>

#include "softplus.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec Logistic(const arma::vec& x) {

  vec res(x.n_elem);

  for (uword i = 0; i < x.n_elem; ++i) {
    res(i) = 1 / (expm1(-x(i)) + 2);
  }

  return res;
}


// [[Rcpp::export]]
arma::vec Softplus1(const arma::vec& x) {

  vec res(x.n_elem);

  for (uword i = 0; i < x.n_elem; ++i) {
    res(i) = log1p(exp(x(i)));
  }

  return res;
}


// [[Rcpp::export]]
arma::vec Softplus(const arma::vec& x,
                   const arma::vec& weight) {

  vec res(x.n_elem);

  for (uword i = 0; i < x.n_elem; ++i) {
    res(i) = log1p(exp(x(i))) * weight(i);
  }

  return res;
}



// [[Rcpp::export]]
arma::vec SoftplusGrad1(const arma::vec& x) {

  return Logistic(x) / sum(Softplus1(x));

}

// [[Rcpp::export]]
arma::vec SoftplusGrad(const arma::vec& x,
                       const arma::vec& weight) {

  return Logistic(x) % weight / sum(Softplus(x, weight));

}


