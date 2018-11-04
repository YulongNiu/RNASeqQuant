#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>

#include "softmax.h"
#include "softplus.h"
#include "gradient.h"
#include "utilities.h"

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

  for (auto i : idx) {
    grad.elem(ec[i]) += count(i) * Softmax(w.elem(ec[i]), 1/efflen[i]);
  }

  // std::cout << grad.subvec(0, 9) << std::endl;

  return sum(count.elem(idx)) * Softmax1(w) - grad;
}


// [[Rcpp::export]]
arma::vec GradientSM2(const arma::vec& w,
                      const std::vector< std::vector< arma::vec > >& efflen,
                      const std::vector< std::vector< arma::uvec > >& ec,
                      const arma::uvec& count,
                      const arma::uvec& spenum,
                      const arma::uvec& idx) {

  // split species number
  uword sn = spenum.n_elem - 1;
  vector<vec> wsplit(sn);
  vec wlse(sn);
  vec wsf(w.n_elem);

  for (uword i = 0; i < sn; ++i) {

    uword start = spenum(i);
    uword end = spenum(i) + spenum(i+1) - 1;
    vec eachw = w.subvec(start, end);
    wsplit[i] = eachw;
    wlse(i) = LogSumExp1(eachw);
    wsf.subvec(start, end) = Softmax1(eachw);
  }

  // split w
  uword ecn = ec.size();
  vector< vector< vec > > ecw(ecn, vector< vec >(sn));
  for (uword i = 0; i < ecn; ++i) {
    for (uword j = 0; j < sn; ++j) {
      ecw[i][j] = w.elem(ec[i][j]);
    }
  }

  // compute gradient
  vec grad(w.n_elem, fill::zeros);
  for (auto i : idx) {
    grad.elem(CmpUvec(ec.at(i))) += count(i) * ECGradSM(wsplit, wlse, efflen[i], ecw[i]);
  }

  return sum(count.elem(idx)) * wsf - grad;
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
    grad.elem(ec[ei]) += count(ei) * SoftplusGrad(w.elem(ec[ei]), 1/efflen[ei]);
  }

  return sum(count.elem(idx)) * SoftplusGrad1(w) - grad;
}
