#include <RcppArmadillo.h>

#include <vector>

#include "utilities.h"
#include "gradient.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec Test(const Rcpp::CharacterVector& ecraw,
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

  vec res = GradientSM2(w, efflen, ec, count, spenum, idx);
  return res;
}
