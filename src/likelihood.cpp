#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double LL(const arma::vec& prob,
          const std::vector<arma::vec>& efflen,
          const std::vector<arma::uvec>& ec,
          const arma::uvec& count) {

  uword ecnum = ec.size();
  vec eachll(ecnum, fill::zeros);

  for (uword i = 0; i < ecnum; ++i) {
    eachll(i) = count(i) * log(sum(prob.elem(ec[i]) / efflen[i]));
  }

  return sum(eachll);
}
