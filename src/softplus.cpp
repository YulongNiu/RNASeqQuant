#include <RcppArmadillo.h>

#include <cmath>

#include "softplus.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


//' Calculate the logistic and softplus calculator
//'
//' \itemize{
//'   \item \code{Logistic()}: The logistic function.
//'   \item \item \code{Softplus()} and \code{Softplus1()}: Softplus with or without weight.
//'   \item \code{SoftplusGrad()} and \code{SoftplusGrad1()}: Internal functions for partial derivation.
//' }
//'
//' @title Softplus
//' @return
//' \itemize{
//'   \item \code{Logisitc()}: A \code{arma::vec} indicates the logistic function \eqn{\frac{1}{1 + \mathrm{e}^{-x}}}.
//'   \item \code{Softplus()} and \code{Softplus1()}: A \code{arma::vec} indicates softplus with (\eqn{\log(1 + \mathrm{e}^{x_i}) * weight_i} or without weight (\eqn{\log(1 + \mathrm{e}^{x_i})}).
//'   \item \code{SoftplusGrad()} and \code{SoftplusGrad1()}: A \code{arma::vec} indicates \eqn{\frac{logistic * weight_i}{\sum{softplus}}} or \eqn{\frac{logistic}{\sum{softplus1}}}.
//' }
//' @inheritParams LogSumExp
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname softplus
//' @keywords internal
// [[Rcpp::export]]
arma::vec Logistic(const arma::vec& x) {

  vec res(x.n_elem);

  for (uword i = 0; i < x.n_elem; ++i) {
    res(i) = 1 / (expm1(-x(i)) + 2);
  }

  return res;
}


//' @inheritParams LogSumExp
//' @rdname softplus
//' @keywords internal
// [[Rcpp::export]]
arma::vec Softplus1(const arma::vec& x) {

  vec res(x.n_elem);

  for (uword i = 0; i < x.n_elem; ++i) {
    res(i) = log1p(exp(x(i)));
  }

  return res;
}


//' @inheritParams LogSumExp
//' @rdname softplus
//' @keywords internal
// [[Rcpp::export]]
arma::vec Softplus(const arma::vec& x,
                   const arma::vec& weight) {

  return Softplus1(x) % weight;

}


//' @inheritParams LogSumExp
//' @rdname softplus
//' @keywords internal
// [[Rcpp::export]]
arma::vec SoftplusGrad1(const arma::vec& x) {

  return Logistic(x) / sum(Softplus1(x));

}


//' @inheritParams LogSumExp
//' @rdname softplus
//' @keywords internal
// [[Rcpp::export]]
arma::vec SoftplusGrad(const arma::vec& x,
                       const arma::vec& weight) {

  return Logistic(x) % weight / sum(Softplus(x, weight));

}


