#include <RcppArmadillo.h>

#include <vector>

using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec InvSqrtRoot(const arma::vec& x,
                      const double alpha) {

  return 1 / sqrt(1 + alpha * x % x);

}

// [[Rcpp::export]]
arma::vec ISRU1(const arma::vec& x,
                const arma::vec& isr,
                const double alpha) {

  return isr % x + 1 / sqrt(alpha);

}

// [[Rcpp::export]]
arma::vec ISRU(const arma::vec& x,
               const arma::vec& isr,
               const arma::vec& weight,
               const double alpha) {

  return  isr % x % weight + weight / sqrt(alpha);

}


// [[Rcpp::export]]
arma::vec ISRUGrad1(const arma::vec& x,
                    const arma::vec& isr,
                    const double alpha) {

  // numerator
  vec nr = pow(isr, 3);

  // denominator
  double dn = sum(isr % x) + x.n_elem / sqrt(alpha);

  return nr / dn;
}


// [[Rcpp::export]]
arma::vec ISRUGrad(const arma::vec& x,
                   const arma::vec& isr,
                   const arma::vec& weight,
                   const double alpha) {

  // numerator
  vec nr = pow(isr, 3) % weight;

  // denominator
  double dn = sum(isr % x % weight) + sum(weight) / sqrt(alpha);

  return nr / dn;
}
