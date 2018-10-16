#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>

#include "utilities.h"
#include "likelihood.h"

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


// [[Rcpp::export]]
arma::vec Estw2Estcount(const arma::vec& estw,
                        double cn) {
  return exp(estw) / sum(exp(estw)) * cn;
}


// [[Rcpp::export]]
arma::vec BGD(const arma::vec& efflenraw,
              const Rcpp::CharacterVector& ecraw,
              const arma::uvec& countraw,
              const arma::uvec& spenumraw,
              const arma::uword maxiter = 10000,
              const arma::uword miniter = 50,
              const double alpha = 0.01) {

  // stop iteration params from kallisto
  double countChangeLimit = 1e-2;
  double countChange = 1e-2;
  double countLimit = 1e-8;

  // step1: pseudo information
  // remove zero counts
  uvec zeros = find(countraw > 0);
  IntegerVector zerosidx(zeros.begin(), zeros.end());

  uvec count = countraw.elem(zeros);
  vector<uvec> ec = SplitEC(ecraw[zerosidx]);
  vector<vec> efflen = MatchEfflen(ec, efflenraw);
  uvec spenum = IdxSpenum(spenumraw);

  // step2: EM iteration
  // start w and estcount
  uword tn = sum(spenumraw);
  double cn = sum(count);
  vec startw(tn, fill::ones);
  vec w(tn, fill::zeros);
  vec startest = Estw2Estcount(w, cn);
  vec est(tn, fill::zeros);
  vec prob(tn, fill::zeros);

  for (uword iter = 0; iter < maxiter; ++iter) {

    w = startw - alpha * Gradient(startw, efflen, ec, count);
    est = Estw2Estcount(w, cn);
    prob = exp(w) / sum(exp(w));

    Rcout << std::setprecision (20) << LL(prob, efflen, ec, count) << std::endl;
    // cout << std::setprecision (20) << sum(est) << endl;
    // cout << sum(prob) << endl;

    // stop iteration condition
    uword nopassn = 0;
    for (uword t = 0; t < tn; ++t) {
      if (est(t) > countChangeLimit && (fabs(est(t) - startest(t))/est(t)) > countChange) {
        ++nopassn;
      } else {}
    }

    if (nopassn == 0 && iter >= miniter) {
      Rcout << "The iteration number is " << iter + 1
            << ". The log likelihood is " << LL(prob, efflen, ec, count)
            << "." << std::endl;
      break;
    } else {
      startw = w;
      startest = est;
    }
  }

  // reset small est
  est.elem(find(est < countLimit)).zeros();

  return est;
}
