#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

// [[Rcpp::export]]
arma::vec Gradient(const arma::vec& w,
                   const std::vector<arma::vec>& efflen,
                   const std::vector<arma::uvec>& ec,
                   const arma::uvec& count) {

  vec grad(w.n_elem, fill::zeros);

  // exponent
  vec ew = exp(w);
  double Z = sum(ew);

  for (uword i = 0; i < count.n_elem; ++i) {
    vec eachewp = ew.elem(ec[i]) / efflen[i];
    grad.elem(ec[i]) += eachewp * count(i) / sum(eachewp);
  }

  return sum(count) * ew / Z - grad;
}
