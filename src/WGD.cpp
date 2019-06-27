#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>

#include "utilities.h"
#include "likelihood.h"
#include "AFfactory.h"
#include "AFmeasure.h"
#include "activation.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]


// [[Rcpp::export]]
arma::uvec CountRepeat(arma::uvec x) {

  vec u(1);
  uvec uNum(1);
  uword xL = x.n_elem;
  uword j;

  u(0) = x(0);
  uNum(0) = 1;

  for(uword i = 1; i < xL; ++i) {

    uword n = u.n_elem;

    for (j = 0; j < n; ++j) {
      if (x(i) == u(j)) {
        uNum(j) += 1;
        break;
      } else {}
    }

    if (j == n) {
      // add new unique elements
      u.resize(n + 1);
      u(n) = x(i);
      uNum.resize(n + 1);
      uNum(n) = 1;
    } else {}

  }

  return uNum;
}


// [[Rcpp::export]]
arma::uvec CountEC(const std::vector<arma::uvec>& ec) {

  // collapse ec
  uvec coles;
  for (auto i : ec) {
    coles = join_cols(coles, i);
  }

  return CountRepeat(coles);
}


// [[Rcpp::export]]
arma::vec AdamW(const arma::vec& efflenraw,
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
  double beta1 = 0.9;
  double beta2 = 0.999;
  double epsilon = 1e-8;

  // step1: pseudo information remove zero counts
  uvec zeros = find(countraw > 0);
  IntegerVector zerosidx(zeros.begin(), zeros.end());

  uvec count = countraw.elem(zeros);
  vector<uvec> ec = SplitEC(ecraw[zerosidx]);
  vector<vec> efflen = MatchEfflen(ec, efflenraw);

  // step2: Adam
  uword tn = sum(spenumraw);
  uword cn = sum(count);
  // uword sn = spenumraw.n_elem;
  uword ecn = ec.size();

  // Glorot normal initializer/Xavier normal initializer
  vec w = randn<vec>(tn) / sqrt(tn);
  // vec w(tn); w.fill(0.01);

  uvec ecw = 1/CountEC(SplitEC(ecraw));
  vec m = vec(tn, fill::zeros);
  vec v = vec(tn, fill::zeros);
  uword t = 0;

  // gradient and shuffled index
  vec grad = vec(tn);
  uvec idx = linspace<uvec>(0, ecn - 1, ecn);

  // active function
  std::shared_ptr<AFmeasure> afgrad = AFfactory().createAFGradient(attrs, arguments);
  std::shared_ptr<AFmeasure> afc = AFfactory().createAFCounts(attrs, arguments);

  for (uword iter = 0; iter < epochs; ++iter) {

    // Rcout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(afc->AFCounts(w), efflen, ec, count) << std::endl;

    idx = shuffle(idx);
    uword biter = 0;

    // mini-batch
    while (biter < ecn) {
      ++t;
      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      // adam for each batch
      grad = afgrad->AFGradient(w, efflen, ec, count, eachidx);
      m = beta1 * m + (1 - beta1) * grad;
      v = beta2 * v + (1 - beta2) * square(grad);
      double etat = eta * sqrt(1 - pow(beta2, t)) / (1 - pow(beta1, t));
      w = w - etat * m / (sqrt(v) + epsilon) - eta * ecw % grad;

      biter += batchsize;
    }
  }

  // step3: reset small est
  vec est = afc->AFCounts(w) * cn;
  Rcout << "The log likelihood is " << std::setprecision (20) << LL(est, efflen, ec, count) <<
    "." << std::endl;

  est.elem(find(est < countLimit)).zeros();

  return est;
}




// [[Rcpp::export]]
arma::vec NRMSPropW(const arma::vec& efflenraw,
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

  // adagrad settings
  double gamma = 0.9;
  double epsilon = 1e-8;
  double velocity = 0.9;

  // step1: pseudo information
  // remove zero counts
  uvec zeros = find(countraw > 0);
  IntegerVector zerosidx(zeros.begin(), zeros.end());

  uvec count = countraw.elem(zeros);
  vector<uvec> ec = SplitEC(ecraw[zerosidx]);
  vector<vec> efflen = MatchEfflen(ec, efflenraw);

  // step2: Adagrad
  // start w and estcount
  uword tn = sum(spenumraw);
  uword cn = sum(count);
  // uword sn = spenumraw.n_elem;
  uword ecn = ec.size();

  // Glorot normal initializer/Xavier normal initializer
  vec w = randn<vec>(tn) / sqrt(tn);
  // vec w(tn); w.fill(0.01);

  uvec ecw = CountEC(SplitEC(ecraw));
  vec eg2 = vec(tn, fill::zeros);
  vec V = vec(tn, fill::zeros);

  // gradient and shuffled index
  vec grad = vec(tn);
  uvec idx = linspace<uvec>(0, ecn - 1, ecn);

  // active function
  std::shared_ptr<AFmeasure> afgrad = AFfactory().createAFGradient(attrs, arguments);
  std::shared_ptr<AFmeasure> afc = AFfactory().createAFCounts(attrs, arguments);

  for (uword iter = 0; iter < epochs; ++iter) {

    // std::cout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(afc->AFCounts(w), efflen, ec, count) << "|" << t << std::endl;
    idx = shuffle(idx);
    uword biter = 0;

    // mini-batch
    while (biter < ecn) {
      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      // NAG
      grad = afgrad->AFGradient(w - velocity * V, efflen, ec, count, eachidx);
      eg2 = gamma * eg2 + (1 - gamma) * grad % grad;

      // update V
      V = velocity * V + eta / sqrt(eg2 + epsilon) % grad + eta * ecw % grad;
      w -= V;

      biter += batchsize;
    }
  }

  // reset small est
  vec est = afc->AFCounts(w) * cn;
  Rcout << "The log likelihood is " << std::setprecision (20) << LL(est, efflen, ec, count) <<
    "." << std::endl;

  est.elem(find(est < countLimit)).zeros();

  return est;
}
