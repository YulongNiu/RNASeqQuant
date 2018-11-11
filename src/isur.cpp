#include <RcppArmadillo.h>

#include <vector>

#include "utilities.h"

using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec SquareRoot(const arma::vec& x,
                     const double alpha) {

  return sqrt(1 + alpha * x * x);

}

// [[Rcpp::export]]
arma::vec ISRU1(const arma::vec& x,
                const double alpha) {

  return x % SquareRoot(x, alpha) + 1 / sqrt(alpha);

}

// [[Rcpp::export]]
arma::vec ISRU(const arma::vec& x,
                const arma::vec& weight,
                const double alpha) {

  return ISRU1(x, alpha) % weight;

}


// [[Rcpp::export]]
arma::vec ISRUGrad1(const arma::vec& x,
                    const double alpha) {

  // numerator
  vec sr = SquareRoot(x, alpha);
  vec nr = pow(sr, 3);

  // denominator
  double dn = sum(sr % x) + x.n_elem / sqrt(alpha);

  return nr / dn;
}


// [[Rcpp::export]]
arma::vec ISRUGrad(const arma::vec& x,
                   const arma::vec& weight,
                   const double alpha) {

  // numerator
  vec sr = SquareRoot(x, alpha);
  vec nr = pow(sr, 3) % weight;

  // denominator
  double dn = sum(sr % x % weight) + sum(weight) / sqrt(alpha);

  return nr / dn;
}



