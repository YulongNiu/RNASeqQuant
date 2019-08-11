#include <RcppArmadillo.h>

#include <algorithm>
#include <vector>

#include "utilities.h"
#include "likelihood.h"
#include "AFfactory.h"
#include "AFmeasure.h"
#include "activation.h"
#include "Optimizer.h"
#include "Optfactory.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
void Test() {

  Optimizer *testadam  = new NRMSProp(10, 0.9, 0.9, 1e-8);

  vec g = vec(10); g.fill(0.1);
  vec w = vec(10); w.fill(20);

  Rcout << testadam -> update(w, g, 0.1) << endl;
  Rcout << testadam -> preupdate(w) << endl;

}


// TestObj(10, list(opt = 'Adam'), list())
// [[Rcpp::export]]
void TestObj(arma::uword tn,
             const Rcpp::List attrs,
             const Rcpp::List arguments) {

  std::shared_ptr<Optimizer> optobj = Optfactory().createOpt(tn, attrs, arguments);

  vec g = vec(10); g.fill(0.1);
  vec w = vec(10); w.fill(20);

  Rcout << optobj -> update(w, g, 0.1) << endl;
}


// [[Rcpp::export]]
arma::vec GD(const arma::vec& efflenraw,
             const Rcpp::CharacterVector& ecraw,
             const arma::uvec& countraw,
             const arma::uvec& spenumraw,
             const arma::uword epochs,
             const arma::uword batchsize,
             const double eta,
             const Rcpp::List attrs,
             const Rcpp::List arguments) {

  // stop iteration settings from kallisto
  // double countChangeLimit = 1e-2
  // double countChange = 1e-2
  double countLimit = 1e-8;

  // adam settings
  double decay = 0.03;

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

  // Glorot normal initializer/Xavier normal initializer
  vec w = randn<vec>(tn) / sqrt(tn);
  // vec w(tn); w.fill(0.01);
  vec grad = vec(tn);
  uvec idx = linspace<uvec>(0, ecn - 1, ecn);

  // active function
  std::shared_ptr<AFmeasure> afgrad = AFfactory().createAFGradient(attrs, arguments);
  std::shared_ptr<AFmeasure> afc = AFfactory().createAFCounts(attrs, arguments);

  // GD method
  std::shared_ptr<Optimizer> gd = Optfactory().createOpt(tn, attrs, arguments);

  // step3: gradient decent
  for (uword iter = 0; iter < epochs; ++iter) {

    // Rcout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(afc->AFCounts(w), efflen, ec, count) << std::endl;

    idx = shuffle(idx);
    uword biter = 0;
    double etai = eta / (1 + decay * iter);

    // mini-batch
    while (biter < ecn) {

      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      grad = afgrad -> AFGradient(gd -> preupdate(w), efflen, ec, count, eachidx);
      w = gd -> update(w, grad, etai);
      biter += batchsize;
    }
  }

  // step3: reset small est
  vec est = afc -> AFCounts(w) * cn;
  Rcout << "The log likelihood is " << std::setprecision (20) << LL(est, efflen, ec, count) <<
    "." << std::endl;

  est.elem(find(est < countLimit)).zeros();

  return est;
}
