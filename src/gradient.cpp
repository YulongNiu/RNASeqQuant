#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>

#include "softmax.h"
#include "softplus.h"
#include "gradient.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

// [[Rcpp::export]]
arma::vec GradientSM(const arma::vec& w,
                     const std::vector<arma::vec>& efflen,
                     const std::vector<arma::uvec>& ec,
                     const arma::uvec& count,
                     const arma::uvec& idx) {

  vec grad(w.n_elem, fill::zeros);

  for (uword i = 0; i < idx.n_elem; ++i) {
    uword ei = idx(i);
    grad.elem(ec[ei]) += count(ei) * Softmax(w.elem(ec[ei]), 1/efflen[ei]);
  }

  return sum(count.elem(idx)) * Softmax1(w) - grad;
}


// [[Rcpp::export]]
arma::vec GradientSP(const arma::vec& w,
                     const std::vector<arma::vec>& efflen,
                     const std::vector<arma::uvec>& ec,
                     const arma::uvec& count,
                     const arma::uvec& idx) {

  vec grad(w.n_elem, fill::zeros);

  for (uword i = 0; i < idx.n_elem; ++i) {
    uword ei = idx(i);
    grad.elem(ec[ei]) += count(ei) * Softplus(w.elem(ec[ei]), 1/efflen[ei]);
  }

  return sum(count.elem(idx)) * Softplus1(w) - grad;
}
