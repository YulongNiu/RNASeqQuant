#include <RcppArmadillo.h>

#include <algorithm>
#include <vector>

#include "utilities.h"
#include "softmax.h"
#include "softplus.h"
#include "isru.h"
#include "gradient.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


//' Gradient.
//'
//' \itemize{
//'   \item \code{GradientSM_()}: Gradients of Softmax.
//'   \item \code{GradientSP_()}: Gradients of SoftPlus.
//'   \item \code{GradientISRU_()}: Gradients of ISRU.
//' }
//'
//' @title Calculate gradient
//' @return A \code{arma::vec} indicates gradients.
//' @param w A \code{arma::vec} indicates estimated parameters.
//' @param idx A \code{arma::uvec} indicates the index of \code{w} used for gradient descending.
//' @inheritParams LL
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname gradient
//' @keywords internal
// [[Rcpp::export]]
arma::vec GradientSM_(const arma::vec& w,
                      const std::vector<arma::vec>& efflen,
                      const std::vector<arma::uvec>& ec,
                      const arma::uvec& count,
                      const arma::uvec& idx) {

  vec grad(w.n_elem, fill::zeros);

  for (auto i : idx) {
    grad.elem(ec[i]) += count(i) * Softmax(w.elem(ec[i]), 1/efflen[i]);
  }

  return sum(count.elem(idx)) * Softmax1(w) - grad;
}

//' @inheritParams GradientSM_
//' @rdname gradient
//' @keywords internal
// [[Rcpp::export]]
arma::vec GradientSP_(const arma::vec& w,
                      const std::vector<arma::vec>& efflen,
                      const std::vector<arma::uvec>& ec,
                      const arma::uvec& count,
                      const arma::uvec& idx) {

  vec grad(w.n_elem, fill::zeros);

  for (auto i : idx) {
    grad.elem(ec[i]) += count(i) * SoftplusGrad(w.elem(ec[i]), 1/efflen[i]);
  }

  return sum(count.elem(idx)) * SoftplusGrad1(w) - grad;
}


//' @inheritParams GradientSM_
//' @inheritParams InvSqrtRoot
//' @rdname gradient
//' @keywords internal
// [[Rcpp::export]]
arma::vec GradientISRU_(const arma::vec& w,
                        const std::vector<arma::vec>& efflen,
                        const std::vector<arma::uvec>& ec,
                        const arma::uvec& count,
                        const arma::uvec& idx,
                        const double alpha) {

  vec grad(w.n_elem, fill::zeros);

  for (auto i : idx) {
    grad.elem(ec[i]) += count(i) * ISRUGrad(w.elem(ec[i]), InvSqrtRoot(w.elem(ec[i]), alpha), 1/efflen[i], alpha);
  }

  return sum(count.elem(idx)) * ISRUGrad1(w, InvSqrtRoot(w, alpha), alpha) - grad;
}


