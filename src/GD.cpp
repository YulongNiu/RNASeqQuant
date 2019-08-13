#include <RcppArmadillo.h>

#include <algorithm>
#include <vector>

#include "utilities.h"
#include "likelihood.h"
#include "AFfactory.h"
#include "Optfactory.h"

using namespace Rcpp;
using namespace arma;
using namespace std;



// [[Rcpp::export]]
Rcpp::List GD(const arma::vec& efflenraw,
              const Rcpp::CharacterVector& ecraw,
              const arma::uvec& countraw,
              const arma::uvec& spenumraw,
              const arma::uword epochs,
              const arma::uword batchsize,
              const Rcpp::List attrs,
              const Rcpp::List arguments,
              const bool details = false) {

  // stop iteration settings from kallisto
  // double countChangeLimit = 1e-2
  // double countChange = 1e-2
  double countLimit = 1e-8;

  // gd settings
  double eta = arguments.containsElementNamed("eta") ? arguments["eta"] : 0.1;
  double decay = arguments.containsElementNamed("decay") ? arguments["decay"] : 0.03;

  // step1: pseudo information remove zero counts
  uvec zeros = find(countraw > 0);
  IntegerVector zerosidx(zeros.begin(), zeros.end());

  uvec count = countraw.elem(zeros);
  vector<uvec> ec = SplitEC(ecraw[zerosidx]);
  vector<vec> efflen = MatchEfflen(ec, efflenraw);

  // step2: initialization
  uword tn = sum(spenumraw);
  uword cn = sum(count);
  uword ecn = ec.size();
  uvec ftidx = FalseTIdx(ec, spenumraw);
  vec resll(epochs, fill::zeros);

  // Glorot normal initializer/Xavier normal initializer
  vec w = randn<vec>(tn) / sqrt(tn);
  // vec w(tn); w.fill(0.01);
  w.elem(ftidx).fill(-1e8);
  vec grad = vec(tn);
  uvec idx = linspace<uvec>(0, ecn - 1, ecn);

  // active function
  std::shared_ptr<AFmeasure> af = AFfactory().createAF(attrs, arguments);

  // GD method
  std::shared_ptr<Optimizer> gd = Optfactory().createOpt(tn, attrs, arguments);

  // step3: gradient decent
  uword iter;
  for (iter = 0; iter < epochs; ++iter) {

    // Rcout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(af->AFCounts(w) * cn, efflen, ec, count) << std::endl;

    if (details) {
      resll(iter) = LL(af->AFCounts(w) * cn, efflen, ec, count);
    } else {}

    idx = shuffle(idx);
    uword biter = 0;
    double etai = eta / (1 + decay * iter);

    // mini-batch
    while (biter < ecn) {

      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      grad = af->AFGradient(gd->preupdate(w), efflen, ec, count, eachidx);
      w = gd->update(w, grad, etai);

      biter += batchsize;
    }
  }

  // step4: small est & no ec transcripts --> zero
  vec est = af->AFCounts(w) * cn;
  est.elem(find(est < countLimit)).zeros();

  List res = List::create(_["counts"] = est,
                          _["ll"] = resll.subvec(0, iter - 1));

  Rcout << "The iteration number is " << epochs
        << ". The log likelihood is " << std::setprecision (20) << LL(est, efflen, ec, count)
        << "." << std::endl;

  return res;
}
