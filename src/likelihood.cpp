#include <RcppArmadillo.h>

#include "utilities.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


//' Logistic likelihood.
//'
//' The log-likelihood of given probabilities.
//'
//' @title Calculate log-likelihood
//' @return A \code{double} indicates log-likelihood.
//' @param prob A \code{arma::vec} indicates probabilities of selecting a read from the different transcripts.
//' @param efflen A \code{std::vector<arma::vec>} indicated effective length of transcripts.
//' @param ec A \code{std::vector<arma::uvec>} indicated equivalence classes (ec).
//' @param count A \code{arma::uvec} indicated counts of ec.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
double LLEM(const arma::vec& prob,
            const std::vector<arma::vec>& efflen,
            const std::vector<arma::uvec>& ec,
            const arma::uvec& count) {

  uword ecnum = ec.size();
  vec eachll(ecnum);

  for (uword i = 0; i < ecnum; ++i) {
    eachll(i) = count(i) * log(sum(prob.elem(ec[i]) / efflen[i]));
  }

  // std::cout << std::setprecision (20) << sum(eachll) << std::endl;

  return sum(eachll);
}


// [[Rcpp::export]]
double LLGD(const arma::vec& prob,
            const std::vector< std::vector< arma::vec > >& efflen,
            const std::vector< std::vector< arma::uvec > >& ec,
            const arma::uvec& count) {

  uword ecnum = ec.size();
  vec eachll(ecnum);

  for (uword i = 0; i < ecnum; ++i) {
    eachll(i) = count(i) * log(sum(prob.elem(CmpUvec(ec[i])) / CmpVec(efflen[i])));
  }

  return sum(eachll);
}
