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

  // std::cout << grad.subvec(0, 9) << std::endl;

  return sum(count.elem(idx)) * Softmax1(w) - grad;
}


// [[Rcpp::export]]
arma::vec GradientSM2(const arma::vec& w,
                      const std::vector<arma::vec>& efflen,
                      const std::vector<arma::uvec>& ec,
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

  // for (auto x : wsplit) {
  //   std::cout << x.subvec(0, 9) << std::endl;
  // }
  // std::cout << wlse << std::endl;
  // std::cout << wsf.subvec(0, 9) << std::endl;

  // compute gradient
  vec grad(w.n_elem, fill::zeros);

  for (uword i = 0; i < idx.n_elem; ++i) {
    uword ei = idx(i);
    grad.elem(ec[ei]) += count(ei) * ECGradSM(wsplit, wlse, efflen[ei], ec[ei], spenum);
  }

  // std::cout << grad.subvec(0, 9) << std::endl;

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
