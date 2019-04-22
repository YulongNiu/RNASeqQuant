#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>

#include "utilities.h"
#include "likelihood.h"
#include "softmax.h"
#include "softplus.h"
#include "isru.h"
#include "gradient.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]


// [[Rcpp::export]]
arma::vec Adam(const arma::vec& efflenraw,
               const Rcpp::CharacterVector& ecraw,
               const arma::uvec& countraw,
               const arma::uvec& spenumraw,
               const arma::uword epochs = 300,
               const arma::uword batchsize = 1000,
               const double alpha = 0.1) {

  // stop iteration settings from kallisto
  // double countChangeLimit = 1e-2
  // double countChange = 1e-2
  double countLimit = 1e-8;

  // adam settings
  double beta1 = 0.9;
  double beta2 = 0.999;
  double epsilon = 1e-8;

  // step1: pseudo information
  // remove zero counts
  uvec zeros = find(countraw > 0);
  IntegerVector zerosidx(zeros.begin(), zeros.end());

  uvec count = countraw.elem(zeros);
  vector<uvec> ec = SplitEC(ecraw[zerosidx]);
  vector<vec> efflen = MatchEfflen(ec, efflenraw);

  // step2: Adam
  // start w and estcount
  uword tn = sum(spenumraw);
  uword cn = sum(count);
  // uword sn = spenumraw.n_elem;
  uword ecn = ec.size();

 // Glorot normal initializer/Xavier normal initializer
  vec w = randn<vec>(tn) / sqrt(tn);
  // vec w(tn); w.fill(0.01);
  vec m = vec(tn, fill::zeros);
  vec v = vec(tn, fill::zeros);
  uword t = 0;
  // gradient and shuffled index
  vec grad = vec(tn);
  uvec idx = linspace<uvec>(0, ecn - 1, ecn);

  for (uword iter = 0; iter < epochs; ++iter) {

    // std::cout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(Softmax1(w), efflen, ec, count) << "|" << t << std::endl;

    // std::cout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(Softplus1(w) / sum(Softplus1(w)), efflen, ec, count) << "|" << t << std::endl;

    // std::cout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(ISRU1(w, InvSqrtRoot(w, alpha), alpha) / sum(ISRU1(w, InvSqrtRoot(w, alpha), alpha)), efflen, ec, count) << "|" << t << std::endl;

    idx = shuffle(idx);
    uword biter = 0;

    // mini-batch
    while (biter < ecn) {
      ++t;
      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      // adam for each batch
      // grad = GradientSM(w, efflen, ec, count, eachidx);
      // grad = GradientSP(w, efflen, ec, count, eachidx);
      grad = GradientISRU(w, efflen, ec, count, alpha, eachidx);
      m = beta1 * m + (1 - beta1) * grad;
      v = beta2 * v + (1 - beta2) * square(grad);
      double alphat = alpha * sqrt(1 - pow(beta2, t)) / (1 - pow(beta1, t));
      w -= alphat * m / (sqrt(v) + epsilon);

      biter += batchsize;
    }
  }


  // Rcout << "The log likelihood is " << std::setprecision (20) << LL(Softmax1(w), efflen, ec, count) << "." << std::endl;
  // Rcout << "The log likelihood is " << std::setprecision (20) << LL(Softplus1(w) / sum(Softplus1(w)), efflen, ec, count) << "." << std::endl;
  Rcout << "The log likelihood is " << std::setprecision (20) << LL(ISRU1(w, InvSqrtRoot(w, alpha), alpha) / sum(ISRU1(w, InvSqrtRoot(w, alpha), alpha)), efflen, ec, count) << "." << std::endl;


  // reset small est
  // vec est = Softmax1(w) / sum(Softmax1(w)) * cn;
  // vec est = Softplus1(w) / sum(Softplus1(w)) * cn;
  vec est = ISRU1(w, InvSqrtRoot(w, alpha), alpha) / sum(ISRU1(w, InvSqrtRoot(w, alpha), alpha)) * cn;
  est.elem(find(est < countLimit)).zeros();

  return est;
}

