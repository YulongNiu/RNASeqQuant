#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <tthread/tinythread.h>
#include <algorithm>
#include <vector>
#include <cmath>

#include "utilities.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

#if RCPP_PARALLEL_USE_TBB

tbb::mutex countMutex;

struct ExpectEC : public Worker
{
  const arma::vec& prob;
  const std::vector<arma::vec>& efflen;
  const std::vector<arma::uvec>& ec;
  const arma::uvec& count;
  arma::vec& estcount;

  ExpectEC(const arma::vec& prob,
           const std::vector<arma::vec>& efflen,
           const std::vector<arma::uvec>& ec,
           const arma::uvec& count,
           arma::vec& estcount)
    : prob(prob), efflen(efflen), ec(ec), count(count), estcount(estcount) {}

  void operator()(std::size_t begin, std::size_t end) {
    countMutex.lock();
    for (std::size_t i = begin; i < end; ++i) {
      vec eachcp = prob.elem(ec[i]) / efflen[i];
      estcount.elem(ec[i]) += eachcp * count(i)/ sum(eachcp);
    }
    countMutex.unlock();
  }
};

arma::vec EMSingle(const arma::vec& prob,
                   const std::vector<arma::vec>& efflen,
                   const std::vector<arma::uvec>& ec,
                   const arma::uvec& count) {

  vec estcount(prob.n_elem, fill::zeros);

  // create parallel worker and call
  ExpectEC expectEC(prob, efflen, ec, count, estcount);

  // grain size is 1k
  parallelFor(0, efflen.size(), expectEC, 1000);

  return estcount;
}

#else

arma::vec EMSingle(const arma::vec& prob,
                   const std::vector<arma::vec>& efflen,
                   const std::vector<arma::uvec>& ec,
                   const arma::uvec& count) {

  vec estcount(prob.n_elem, fill::zeros);

  for (uword i = 0; i < count.n_elem; ++i) {
    vec eachcp = prob.elem(ec[i]) / efflen[i];
    estcount.elem(ec[i]) += eachcp * count(i)/ sum(eachcp);
  }

  return estcount;
}

#endif

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


//' Expectation maximization (EM) model for RNA-seq quantification.
//'
//' EM model for RNA-seq quantification. The equivalence class (ec) with 0 counts are removed, because these counts have no contributes to the final results.
//'
//' @title EM model
//' @param countraw A \code{arma::uvec} indicates the counts of ec.
//' @param maxiter The maximum iteration number with the default value of 10000.
//' @param miniter The minimum iteration number with the default value of 50.
//' @inheritParams MatchEfflen
//' @inheritParams SplitEC
//' @inheritParams IdxSpenum
//' @references \href{https://arxiv.org/abs/1104.3889}{Lior Pachter: Models for transcript quantification from RNA-Seq}
//' @return A \code{numeric vector} indicates probabilities of selecting a read from the different transcripts.
//' @examples
//' ## Single species
//' ##    f1 f2 f3
//' ## ec1 1 1 1
//' ## ec2 0 1 1
//' ## ec3 1 0 1
//' ## ec4 1 0 0
//' ## ec5 1 1 0
//' plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))
//' EM(plist$efflen, plist$ec, plist$count, spenum = 3)
//'
//' ## Two species
//' ##    f1 f2 f3 f1' f2'
//' ## ec1 1  1  0  0  1
//' ## ec2 1  0  1  1  0
//' ## ec3 0  1  1  0  0
//' ## ec4 0  0  0  1  1
//' ## ec5 1  0  1  0  1
//' ## ec6 1  1  0  0  0
//' plist <- list(ec = c('0,1,4', '0,2,3', '1,2', '3,4', '0,2,4', '0,1'),
//'               count = rep(1, 6), efflen = rep(1, 5))
//' EM(plist$efflen, plist$ec, plist$count, c(3, 2))
//' ## compare with single species
//' EM(plist$efflen, plist$ec, plist$count, 5)
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @export
// [[Rcpp::export]]
arma::vec EM(const arma::vec& efflenraw,
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

    est = EMSingle(prob, efflen, ec, count);

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
