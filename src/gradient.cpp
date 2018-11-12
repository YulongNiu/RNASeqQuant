#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>

#include "utilities.h"
#include "softmax.h"
#include "softplus.h"
#include "isru.h"
#include "gradient.h"


using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]


// [[Rcpp::export]]
arma::vec GradientSM(const arma::vec& w,
                      const std::vector< std::vector< arma::vec > >& efflen,
                      const std::vector< std::vector< arma::uvec > >& ec,
                      const arma::uvec& count,
                      const arma::uvec& spenum,
                      const arma::uvec& idx) {

  // split species number
  uword sn = spenum.n_elem - 1;
  vector<vec> wsplit(sn);
  vec wlse(sn);
  vec wsm(w.n_elem);

  for (uword i = 0; i < sn; ++i) {
    uword start = spenum(i);
    uword end = spenum(i) + spenum(i+1) - 1;
    vec eachw = w.subvec(start, end);
    wsplit[i] = eachw;
    wlse(i) = LogSumExp1(eachw);
    wsm.subvec(start, end) = Softmax1(eachw);
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

  return sum(count.elem(idx)) * wsm - grad;
}


// [[Rcpp::export]]
arma::vec GradientISRU(const arma::vec& w,
                       const std::vector< std::vector< arma::vec > >& efflen,
                       const std::vector< std::vector< arma::uvec > >& ec,
                       const arma::uvec& count,
                       const arma::uvec& spenum,
                       const double alpha,
                       const arma::uvec& idx) {

  // split species number
  uword sn = spenum.n_elem - 1;
  vector<vec> wsplit(sn);
  vec wlse(sn);
  vec wsm(w.n_elem);

  for (uword i = 0; i < sn; ++i) {
    uword start = spenum(i);
    uword end = spenum(i) + spenum(i+1) - 1;
    vec eachw = w.subvec(start, end);
    wsplit[i] = eachw;
    vec eachisr = InvSqrtRoot(eachw, alpha);
    wlse(i) = ISRU1(eachw, eachisr, alpha);
    wsm.subvec(start, end) = ISRUGrad1(eachw, eachisr, alpha);
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
    grad.elem(CmpUvec(ec.at(i))) += count(i) * ECGradISRU(wsplit, wlse, efflen[i], ecw[i], alpha);
  }

  return sum(count.elem(idx)) * wsm - grad;
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
