#include <RcppArmadillo.h>

#include <vector>

#include "softmax.h"
#include "isru.h"
#include "utilities.h"
#include "gradient.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec GradientSMSS(const arma::vec& w,
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


// [[Rcpp::export]]
arma::vec TestGradientSM(const Rcpp::CharacterVector& ecraw,
                         const arma::vec& efflenraw,
                         const arma::uvec& spenum,
                         const arma::vec& w,
                         const arma::uvec& count,
                         const arma::uvec& idx) {

  vector< vector< uvec > > ec(ecraw.size(), vector< uvec >(spenum.n_elem - 1));
  vector< vector< vec > > efflen(ecraw.size(), vector< vec >(spenum.n_elem - 1));

  EC2Spe(ec, efflen, ecraw, efflenraw, spenum);

  // for (auto x : ec) {
  //   for (auto y : x) {
  //     std::cout << y << std::endl;
  //   }
  // }

  // for (auto x : efflen) {
  //   for (auto y : x) {
  //     std::cout << y << std::endl;
  //   }
  // }

  vec res = GradientSM(w, efflen, ec, count, spenum, idx);
  return res;
}


// [[Rcpp::export]]
arma::vec GradientISRUSS(const arma::vec& w,
                         const std::vector<arma::vec>& efflen,
                         const std::vector<arma::uvec>& ec,
                         const arma::uvec& count,
                         const double alpha,
                         const arma::uvec& idx) {

  vec grad(w.n_elem, fill::zeros);

  for (auto i : idx) {
    grad.elem(ec[i]) += count(i) * ISRUGrad(w.elem(ec[i]), InvSqrtRoot(w.elem(ec[i]), alpha), 1/efflen[i], alpha);
  }

  return sum(count.elem(idx)) * ISRUGrad1(w, InvSqrtRoot(w, alpha), alpha) - grad;
}


// [[Rcpp::export]]
arma::vec TestGradientISRU(const Rcpp::CharacterVector& ecraw,
                           const arma::vec& efflenraw,
                           const arma::uvec& spenum,
                           const arma::vec& w,
                           const arma::uvec& count,
                           const double alpha,
                           const arma::uvec& idx) {

  vector< vector< uvec > > ec(ecraw.size(), vector< uvec >(spenum.n_elem - 1));
  vector< vector< vec > > efflen(ecraw.size(), vector< vec >(spenum.n_elem - 1));

  EC2Spe(ec, efflen, ecraw, efflenraw, spenum);

  vec res = GradientISRU(w, efflen, ec, count, spenum, alpha, idx);
  return res;
}
