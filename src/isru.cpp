#include <RcppArmadillo.h>

#include <vector>

#include "utilities.h"

using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec InvSqrtRoot(const arma::vec& x,
                     const double alpha) {

  return 1/sqrt(1 + alpha * x % x);

}

// [[Rcpp::export]]
double ISRU1(const arma::vec& x,
             const arma::vec& isr,
             const double alpha) {

  return sum(isr % x) + x.n_elem / sqrt(alpha);

}

// [[Rcpp::export]]
double ISRU(const arma::vec& x,
            const arma::vec& isr,
            const arma::vec& weight,
            const double alpha) {

  return  sum(isr % x % weight) + sum(weight) / sqrt(alpha);

}


// [[Rcpp::export]]
arma::vec ISRUGrad1(const arma::vec& x,
                    const double alpha) {

  // numerator
  vec sr = InvSqrtRoot(x, alpha);
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
  vec sr = InvSqrtRoot(x, alpha);
  vec nr = pow(sr, 3) % weight;

  // denominator
  double dn = sum(sr % x % weight) + sum(weight) / sqrt(alpha);

  return nr / dn;
}


// [[Rcpp::export]]
arma::vec SingleSpeGradISRU(const arma::vec& wlse,
                            const std::vector< arma::vec >& efflensg,
                            const std::vector< arma::vec >& ecd,
                            const arma::vec& ecsum,
                            const arma::vec& ecratio,
                            const arma::uword idx) {

  // numerator
  // wratio is 0 if one species has no transcripts in the ec
  // if tratio == 0, then only one species
  double tratio = sum(ecratio) - ecratio(idx);
  vec nr = ecd[idx] % (1 / efflensg[idx] + tratio);

  // denominator
  double dn = ecsum(idx);
  if (tratio > 0) {
    dn += wlse(idx) * tratio;
  } else {}

  return nr / dn;
}


// ECGradISRU(list(c(1, 1), 1), c(2*sqrt(100/101) + 20, sqrt(100/101) + 10), list(1, 1), list(1, 1), 1/100)
// [[Rcpp::export]]
arma::vec ECGradISRU(const std::vector< arma::vec >& w,
                     const arma::vec& wlse,
                     const std::vector< arma::vec >& efflensg,
                     const std::vector< arma::vec >& wsg,
                     const double alpha) {

  // initialization
  uword sn = wsg.size();
  vector<vec> ecd(sn);
  vec ecsum(sn, fill::zeros);
  vec ecratio(sn, fill::zeros);

  // split each ec
  for (uword i = 0; i < sn; ++i) {
    vec eachw = wsg[i];
    if (eachw.n_elem > 0) {
      vec eachisr = InvSqrtRoot(eachw, alpha);
      ecd[i] = pow(eachisr, 3);
      ecsum(i) = ISRU(eachw, eachisr, 1 / efflensg[i], alpha);
      ecratio(i) = ecsum(i) / wlse(i);
    } else {}
  }

  // calculate each species
  vec res;
  for (uword i = 0; i < sn; ++i) {
    if (wsg[i].n_elem > 0) {
      res = join_cols(res, SingleSpeGradISRU(wlse, efflensg, ecd, ecsum, ecratio, i));
    } else {}
  }

  return res;
}


