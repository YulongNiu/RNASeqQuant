#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>

#include "utilities.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

struct ExpectEC : public Worker
{
  arma::vec& prob;
  const std::vector<arma::vec>& efflen;
  const std::vector<arma::uvec>& ec;
  arma::uvec& count;
  arma::vec& estcount;

  ExpectEC(arma::vec& prob,
           const std::vector<arma::vec>& efflen,
           const std::vector<arma::uvec>& ec,
           arma::uvec& count,
           arma::vec& estcount)
    : prob(prob), efflen(efflen), ec(ec), count(count), estcount(estcount) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      vec eachcp = prob.elem(ec[i]) / efflen[i];
      estcount.elem(ec[i]) += eachcp * count(i)/ sum(eachcp);
    }
  }
};


//' Parallel a single EM iteration.
//'
//' One iteration of EM model.
//'
//' @title Single EM iteration
//' @return A updated \code{arma::vec} estimated counts of different transcripts.
//' @param prob A \code{arma::vec} of probabilities of selecting a read from the different transcripts.
//' @param efflen A \code{std::vector<arma::vec>} indicated the effective length of transcripts.
//' @param ec A \code{std::vector<arma::vec>} indicated equivalence classes with the same length of \code{efflen}.
//' @param count A \code{arma::uvec} indicated the mapped count number of equivalence class with the same length of \code{efflen}.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
arma::vec EMSingle(arma::vec& prob,
                   const std::vector<arma::vec>& efflen,
                   const std::vector<arma::uvec>& ec,
                   arma::uvec& count) {

  vec estcount(prob.n_elem, fill::zeros);

  // create parallel worker and call
  ExpectEC expectEC(prob, efflen, ec, count, estcount);
  parallelFor(0, efflen.size(), expectEC);

  return estcount;
}



// [[Rcpp::export]]
arma::vec EMSingle2(arma::vec& prob,
                    const std::vector<arma::vec>& efflen,
                    const std::vector<arma::uvec>& ec,
                    arma::uvec& count) {

  vec estcount(prob.n_elem, fill::zeros);

  for (uword i = 0; i < count.n_elem; ++i) {
    vec eachcp = prob.elem(ec[i]) / efflen[i];
    estcount.elem(ec[i]) += eachcp * count(i)/ sum(eachcp);
  }

  return estcount;
}




//' Transform estimated count to probabilities.
//'
//' Use estimated counts as the outputs and EM stop conditions.
//'
//' @title Counts to probabilities
//' @return A \code{arma::vec} indicates probabilities of selecting a read from the different transcripts.
//' @param estcount A \code{arma::vec} estimated counts of transcripts.
//' @param spenum A \code{arma::uvec} indicated number of transcripts.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
arma::vec Estcount2Prob(const arma::vec& estcount,
                        const arma::uvec& spenum) {

  vec prob(estcount.n_elem, fill::zeros);

  for (uword i = 0; i < spenum.n_elem - 1; ++i) {
    uword start = spenum(i);
    uword end = spenum(i) + spenum(i+1) - 1;
    prob.subvec(start, end) = estcount.subvec(start, end) / sum(estcount.subvec(start, end));
  }

  // // reduce small number influence
  // prob *= 1e+8;

  return prob;
}


// [[Rcpp::export]]
arma::vec EMTest(const arma::vec& efflenraw,
                 const Rcpp::CharacterVector& ecraw,
                 const arma::uvec& countraw,
                 const arma::uvec& spenumraw,
                 const arma::uword maxiter = 10000,
                 const arma::uword miniter = 50) {

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
  // start prob and est
  uword tn = sum(spenumraw);
  double cn = sum(count);
  vec prob(tn);
  vec startest(tn);
  vec est(tn, fill::zeros);
  for (uword i = 0; i < spenum.n_elem - 1; ++i) {
    uword start = spenum(i);
    uword end = spenum(i) + spenum(i+1) - 1;
    prob.subvec(start, end).fill(1.0/spenum(i+1));
    startest.subvec(start, end).fill(cn/(spenum(i+1) * spenumraw.size()));
  }

  for (uword iter = 0; iter < maxiter; ++iter) {

    est = EMSingle2(prob, efflen, ec, count);

    cout << std::setprecision (20) << sum(est) << endl;
    cout << sum(prob) << endl;
    // for (auto eachest : est) {
    //   cout << std::setprecision (20) << eachest << endl;
    // }

    // stop iteration condition
    uword nopassn = 0;
    for (uword t = 0; t < tn; ++t) {
      if (est(t) > countChangeLimit && (fabs(est(t) - startest(t))/est(t)) > countChange) {
        ++nopassn;
        // std::cout << fabs(est(t) - startest(t))/est(t) << endl;
      } else {}
    }

    if (nopassn == 0) {
      std::cout << "The iteration number is " << iter << std::endl;
      break;
    } else {
      prob = Estcount2Prob(est, spenum);
      startest = est;
    }
  }

  // reset small est
  est.elem(find(est < countLimit)).zeros();

  return est;
}

