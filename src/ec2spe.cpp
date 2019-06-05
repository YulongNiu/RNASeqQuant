#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <vector>

#include "ec2spe.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

// [[Rcpp::export]]
arma::uvec AccuSum(const arma::uvec& spenum) {

  uword sn = spenum.n_elem;
  uvec res(sn);

  // use index instead
  res(0) = spenum(0) - 1;

  for (uword i = 0; i < (sn - 1); ++i) {
    res(i + 1) = res(i) + spenum(i + 1);
  }

  return res;
}


// [[Rcpp::export]]
arma::uvec EC2SpeEach(const arma::uvec& accuIdx,
                      const arma::uvec& ec) {

  uword ecn = ec.n_elem;
  uvec eachIdx(ecn);

  for (uword i = 0; i < ecn; ++i) {
    // Rcout << accuIdx >= ec(i) << std::endl;
    eachIdx(i) = find(accuIdx >= ec(i)).min();
  }

  return eachIdx;
}


// [[Rcpp::export]]
std::vector<arma::vec> EC2Spe(const arma::uvec& accuIdx,
                              const std::vector<arma::uvec>& ec,
                              const arma::vec& spefixcounts) {
  uword ecn = ec.size();
  vector<vec> res(ecn);

  for (uword i = 0; i < ecn; ++i) {
    uvec eachidx = EC2SpeEach(accuIdx, ec[i]);
    res[i] = spefixcounts.elem(eachidx);
  }

  return res;
}

