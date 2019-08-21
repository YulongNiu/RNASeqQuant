#include <RcppArmadillo.h>

#include <algorithm>
#include <vector>

#include "utilities.h"
#include "likelihood.h"
#include "AFfactory.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
Rcpp::List RK4(const arma::vec& efflenraw,
                const Rcpp::CharacterVector& ecraw,
                const arma::uvec& countraw,
                const arma::uvec& spenum,
                const Rcpp::List attrs,
                const Rcpp::List arguments,
                const arma::uword maxiter = 10000,
                const arma::uword miniter = 50,
                const arma::uword batchsize = 1024,
                const bool details = false) {

  // stop iteration settings from kallisto
  double countChangeLimit = 1e-2;
  double countChange = 1e-2;
  double countLimit = 1e-8;

  // step1: pseudo information remove zero counts
  uvec zeros = find(countraw > 0);
  IntegerVector zerosidx(zeros.begin(), zeros.end());

  uvec count = countraw.elem(zeros);
  vector<uvec> ec = SplitEC(ecraw[zerosidx]);
  vector<vec> efflen = MatchEfflen(ec, efflenraw);

  // step2: initialization
  uword tn = sum(spenum);
  uword cn = sum(count);
  uword ecn = ec.size();
  vec resll(maxiter, fill::zeros);

  //~~~~~~~~~~~~~~~~~old init~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Glorot normal initializer/Xavier normal initializer
  // vec w(tn); w.fill(0.01);
  // vec w = randn<vec>(tn) / sqrt(tn);
  // uvec ftidx = FalseTIdx(ec, spenum);
  // w.elem(ftidx).fill(-1e8);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  vec w = PreInit(efflen, ec, count, spenum);

  // gd settings
  if (arguments.containsElementNamed("assign0") && !arguments["assign0"]) {
    w.elem(find(w < -1e7)).fill(0.01);
  } else {}
  double eta = arguments.containsElementNamed("eta") ? arguments["eta"] : 0.1;
  double decay = arguments.containsElementNamed("decay") ? arguments["decay"] : 0.03;

  // active function
  std::shared_ptr<AFmeasure> af = AFfactory().createAF(attrs, arguments);
  vec startest = af->AFCounts(w) * cn;
  vec est(tn, fill::zeros);

  // step3: gradient decent
  vec grad = vec(tn);
  uvec idx = linspace<uvec>(0, ecn - 1, ecn);
  uword iter;
  for (iter = 0; iter < maxiter; ++iter) {
    Rcout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(af->AFCounts(w) * cn, efflen, ec, count) << std::endl;

    // record running details
    if (details) {
      resll(iter) = LL(startest, efflen, ec, count);
    } else {}

    idx = shuffle(idx);
    uword biter = 0;
    double etai = eta / (1 + decay * iter);

    // mini-batch
    while (biter < ecn) {

      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      vec k1 = af->AFGradient(w, efflen, ec, count, eachidx);
      vec k2 = af->AFGradient(w + etai/2 * k1, efflen, ec, count, eachidx);
      vec k3 = af->AFGradient(w + etai/2 * k2, efflen, ec, count, eachidx);
      vec k4 = af->AFGradient(w + etai * k3, efflen, ec, count, eachidx);
      w = w - etai * (k1 + 2 * k2 + 2 * k3 + k4) / 6;

      biter += batchsize;
    }

    est = af->AFCounts(w) * cn;

    uword nopassn = sum((est > countChangeLimit) %
                        (abs(est - startest) > countChange));

    // Rcout << nopassn << endl;

    if (nopassn == 0 && iter >= miniter - 1) {
      break;
    } else {
      startest = est;
    }
  }

  // check if maxiter
  if (iter == maxiter) {--iter;} else {}

  // step4: small est & no ec transcripts --> zero
  est.elem(find(est < countLimit)).zeros();

  List res = List::create(_["counts"] = est,
                          _["ll"] = resll.subvec(0, iter));

  Rcout << "The iteration number is " << iter + 1
        << ". The log likelihood is " << std::setprecision (20) << LL(est, efflen, ec, count)
        << "." << std::endl;

  return res;
}
