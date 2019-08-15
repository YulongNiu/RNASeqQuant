#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <tthread/tinythread.h>
#include <algorithm>
#include <vector>
#include <cmath>

#include "utilities.h"
#include "likelihood.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

#if RCPP_PARALLEL_USE_TBB

tbb::mutex countMutex;

struct ExpectEC : public Worker
{
  const arma::vec& startest;
  const std::vector<arma::vec>& efflen;
  const std::vector<arma::uvec>& ec;
  const arma::uvec& count;
  arma::vec& estcount;

  ExpectEC(const arma::vec& startest,
           const std::vector<arma::vec>& efflen,
           const std::vector<arma::uvec>& ec,
           const arma::uvec& count,
           arma::vec& estcount)
    : startest(startest), efflen(efflen), ec(ec), count(count), estcount(estcount) {}

  void operator()(std::size_t begin, std::size_t end) {
    countMutex.lock();
    for (std::size_t i = begin; i < end; ++i) {
      vec eachcp = startest.elem(ec[i]) / efflen[i];
      estcount.elem(ec[i]) += eachcp * count(i)/ sum(eachcp);
    }
    countMutex.unlock();
  }
};

arma::vec EMSingle(const arma::vec& startest,
                   const std::vector<arma::vec>& efflen,
                   const std::vector<arma::uvec>& ec,
                   const arma::uvec& count) {

  vec estcount(startest.n_elem, fill::zeros);

  // create parallel worker and call
  ExpectEC expectEC(startest, efflen, ec, count, estcount);

  // grain size is 1k
  parallelFor(0, efflen.size(), expectEC, 1000);

  return estcount;
}

#else

arma::vec EMSingle(const arma::vec& startest,
                   const std::vector<arma::vec>& efflen,
                   const std::vector<arma::uvec>& ec,
                   const arma::uvec& count) {

  // startest.n_elem is the number of transcripts
  // ec.n_elem == count.n_elem == efflen.n_elem is TRUE, which is the number of equivalent classes/reads
  vec estcount(startest.n_elem, fill::zeros);

  for (uword i = 0; i < ec.size(); ++i) {
    vec eachcp = startest.elem(ec[i]) / efflen[i];
    estcount.elem(ec[i]) += eachcp * count(i)/ sum(eachcp);
  }

  return estcount;
}

#endif


//' Expectation maximization (EM) model for RNA-seq quantification.
//'
//' EM model for RNA-seq quantification. The equivalence class (ec) with 0 counts are removed, because these counts have no contributes to the final results.
//'
//' @title EM model
//' @param countraw A \code{arma::uvec} indicates the counts of ec.
//' @param maxiter The maximum iteration number with the default value of 10000.
//' @param miniter The minimum iteration number with the default value of 50.
//' @param details A \code{bool} value.  When it is set as \code{true}, logistic likelihood and counts for each species in every iteration will be returned, otherwise \code{false}.
//' @inheritParams MatchEfflen
//' @inheritParams SplitEC
//' @inheritParams SpeCount
//' @references \href{https://arxiv.org/abs/1104.3889}{Lior Pachter: Models for transcript quantification from RNA-Seq}
//' @return A \code{List} indicates estimated counts of transcripts.
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
Rcpp::List EM(const arma::vec& efflenraw,
              const Rcpp::CharacterVector& ecraw,
              const arma::uvec& countraw,
              const arma::uvec& spenum,
              const arma::uword maxiter = 10000,
              const arma::uword miniter = 50,
              const bool details = false) {

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

  // step2: EM iteration
  // startest and est
  uword tn = sum(spenum);
  double cn = sum(count);
  uword sn = spenum.n_elem;

  // // cn / tn
  // vec startest(cn/tn);

  // average for each species
  vec scounts(sn);
  vec startest = InitAve(spenum, scounts.fill(cn / sn));
  vec est(tn, fill::zeros);

  // details init
  mat specounts(maxiter, sn, fill::zeros);
  vec resll(maxiter, fill::zeros);

  uword iter;
  for (iter = 0; iter < maxiter; ++iter) {

    // record running details
    if (details) {
      vec eachc = SpeCount(startest, spenum);
      specounts.row(iter) = rowvec(eachc.begin(), sn, false);
      resll(iter) = LL(startest, efflen, ec, count);
    } else {}

    est = EMSingle(startest, efflen, ec, count);

    // stop iteration condition
    uword nopassn = sum((est > countChangeLimit) %
                        ((abs(est - startest) / est) > countChange));

    // Rcout << nopassn << "|" << iter << endl;

    if (nopassn == 0 && iter >= miniter - 1) {
      break;
    } else {
      startest = est;
    }
  }

  // check if maxiter
  if (iter == maxiter) {--iter;} else {}

  // step3: reset small est
  est.elem(find(est < countLimit)).zeros();

  // step4: may add details
  List res = List::create(_["counts"] = est,
                          _["specounts"] = specounts.rows(0, iter),
                          _["ll"] = resll.subvec(0, iter));

  Rcout << "The iteration number is " << iter + 1
        << ". The log likelihood is " << std::setprecision (20) << LL(est, efflen, ec, count)
        << "." << std::endl;

  return res;
}


// [[Rcpp::export]]
Rcpp::List EMSpe(const arma::vec& efflenraw,
                 const Rcpp::CharacterVector& ecraw,
                 const arma::uvec& countraw,
                 const arma::uvec& spenum,
                 const arma::vec& spefixcounts,
                 const arma::uword maxiter = 10000,
                 const arma::uword miniter = 50,
                 const bool details = false) {

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

  // step2: EM iteration
  // startest and est
  uword tn = sum(spenum);
  uword sn = spenum.n_elem;

  // average for each species
  vec startest = InitAve(spenum, spefixcounts);
  vec est(tn, fill::zeros);

  // details init
  mat specounts(maxiter, sn, fill::zeros);
  vec resll(maxiter, fill::zeros);

  uword iter;
  for (iter = 0; iter < maxiter; ++iter) {

    // record running details
    if (details) {
      vec eachc = SpeCount(startest, spenum);
      specounts.row(iter) = rowvec(eachc.begin(), sn, false);
      resll(iter) = LL(startest, efflen, ec, count);
    } else {}

    vec eachlambda = EMSingle(startest, efflen, ec, count);
    est = LambdaSpe(eachlambda, spenum, spefixcounts);

    // stop iteration condition
    uword nopassn = sum((est > countChangeLimit) %
                        ((abs(est - startest) / est) > countChange));

    if (nopassn == 0 && iter >= miniter - 1) {
      Rcout << "The iteration number is " << iter + 1
            << ". The log likelihood is " << std::setprecision (20) << LL(startest, efflen, ec, count)
            << "." << std::endl;
      break;
    } else {
      startest = est;
    }
  }

  // check if maxiter
  if (iter == maxiter) {iter--;} else {}

  // step3: reset small est
  est.elem(find(est < countLimit)).zeros();

  // step4: may add details
  List res = List::create(_["counts"] = est,
                          _["specounts"] = specounts.rows(0, iter),
                          _["ll"] = resll.subvec(0, iter));

  return res;
}

