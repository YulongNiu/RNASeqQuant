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
arma::vec Momentum(const arma::vec& efflenraw,
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
  vec v = vec(tn, fill::zeros);

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
      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      // adam for each batch
      grad = afgrad->AFGradient(w, efflen, ec, count, eachidx);
      v = gamma * v + eta * grad;
      w -= v;

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

// [[Rcpp::export]]
arma::vec NAG(const arma::vec& efflenraw,
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
  vec v = vec(tn, fill::zeros);

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
      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      // adam for each batch
      grad = afgrad->AFGradient(w - gamma * v, efflen, ec, count, eachidx);
      v = gamma * v + eta * grad;
      w -= v;

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


// [[Rcpp::export]]
arma::vec Adam(const arma::vec& efflenraw,
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
  double decay = 0.03;

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
  vec m = vec(tn, fill::zeros);
  vec v = vec(tn, fill::zeros);
  uword t = 0;
  vec grad = vec(tn);
  uvec idx = linspace<uvec>(0, ecn - 1, ecn);

  // active function
  std::shared_ptr<AFmeasure> afgrad = AFfactory().createAFGradient(attrs, arguments);
  std::shared_ptr<AFmeasure> afc = AFfactory().createAFCounts(attrs, arguments);

  for (uword iter = 0; iter < epochs; ++iter) {

    // Rcout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(afc->AFCounts(w), efflen, ec, count) << std::endl;

    idx = shuffle(idx);
    uword biter = 0;
    double etai = eta / (1 + decay * iter);

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
      double etat = etai * sqrt(1 - pow(beta2, t)) / (1 - pow(beta1, t));
      w -= etat * m / (sqrt(v) + epsilon);

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
arma::vec NAdam(const arma::vec& efflenraw,
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
  double decay = 0.03;

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
  vec m = vec(tn, fill::zeros);
  vec v = vec(tn, fill::zeros);
  uword t = 0;
  vec grad = vec(tn);
  uvec idx = linspace<uvec>(0, ecn - 1, ecn);

  // active function
  std::shared_ptr<AFmeasure> afgrad = AFfactory().createAFGradient(attrs, arguments);
  std::shared_ptr<AFmeasure> afc = AFfactory().createAFCounts(attrs, arguments);

  for (uword iter = 0; iter < epochs; ++iter) {

    // std::cout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(afc->AFCounts(w), efflen, ec, count) << "|" << t << std::endl;
    idx = shuffle(idx);
    uword biter = 0;
    double etai = eta / (1 + decay * iter);

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
      vec etat = beta1 * m + (1 - beta1) * grad / (1 - pow(beta1, t));
      w -= etai * etat / (sqrt(v) + epsilon);

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
arma::vec AdaMax(const arma::vec& efflenraw,
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

  vec m = vec(tn, fill::zeros);
  vec u = vec(tn, fill::zeros);
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
      u = Max(beta2 * u, abs(grad));
      w -= eta * m / u;

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
arma::vec Adagrad(const arma::vec& efflenraw,
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
  double epsilon = 1e-8;

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
  vec G = vec(tn, fill::zeros);
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

      // adam for each batch
        grad = afgrad->AFGradient(w, efflen, ec, count, eachidx);
        G += grad % grad;
        w -= eta / sqrt(G + epsilon) % grad;

      biter += batchsize;
    }
  }

  // reset small est
  vec est = afc->AFCounts(w) * cn;
  Rcout << "The log likelihood is " << std::setprecision (20) << LL(afc->AFCounts(w), efflen, ec, count) <<
    "." << std::endl;

  est.elem(find(est < countLimit)).zeros();

  return est;
}


// [[Rcpp::export]]
arma::vec NAdagrad(const arma::vec& efflenraw,
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
  double epsilon = 1e-8;
  double velocity = 0.9;
  double decay = 0.00001;

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
  vec G = vec(tn, fill::zeros);
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
    double etai = eta / (1 + decay * iter);

    // mini-batch
    while (biter < ecn) {
      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      // NAG
      grad = afgrad->AFGradient(w - velocity * V, efflen, ec, count, eachidx);
      G += grad % grad;

      // update V
      V = velocity * V + etai / sqrt(G + epsilon) % grad;
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


// [[Rcpp::export]]
arma::vec Adadelta(const arma::vec& efflenraw,
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
  double gamma = 0.01;
  double epsilon = 1e-8;

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
  vec eg2 = vec(tn, fill::zeros);
  vec edx2 = vec(tn, fill::zeros);

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

      // adam for each batch
      grad = afgrad->AFGradient(w, efflen, ec, count, eachidx);
      eg2 = gamma * eg2 + (1 - gamma) * grad % grad;
      vec dx = -sqrt(edx2 + epsilon) / sqrt(eg2 + epsilon) % grad;
      edx2 = gamma * edx2 + (1 - gamma) * dx % dx;
      w += dx;

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


// [[Rcpp::export]]
arma::vec RMSProp(const arma::vec& efflenraw,
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

  // step1: pseudo information
  // remove zero counts
  uvec zeros = find(countraw > 0);
  IntegerVector zerosidx(zeros.begin(), zeros.end());

  uvec count = countraw.elem(zeros);
  vector<uvec> ec = SplitEC(ecraw[zerosidx]);
  vector<vec> efflen = MatchEfflen(ec, efflenraw);

  // step2: Adagrad
  // start w and estcount
  uvec truetidx = TrueTIdx(ec);
  uword tn = sum(spenumraw);
  uword cn = sum(count);
  // uword sn = spenumraw.n_elem;
  uword ecn = ec.size();

  // Glorot normal initializer/Xavier normal initializer
  vec w = randn<vec>(tn) / sqrt(tn);
  // vec w(tn); w.fill(0.01);
  vec eg2 = vec(tn, fill::zeros);
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

      // adam for each batch
      grad = afgrad->AFGradient(w, efflen, ec, count, eachidx);
      eg2 = gamma * eg2 + (1 - gamma) * grad % grad;
      w -= eta / sqrt(eg2 + epsilon) % grad;

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



// [[Rcpp::export]]
arma::vec NRMSProp(const arma::vec& efflenraw,
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
  double decay = 0.003;

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
  uvec ftidx = FalseTIdx(ec, spenumraw);

  // Glorot normal initializer/Xavier normal initializer
  vec w = randn<vec>(tn) / sqrt(tn);
  // vec w(tn); w.fill(0.01234);
  w.elem(ftidx).fill(-1e8);
  vec eg2 = vec(tn, fill::zeros);
  vec V = vec(tn, fill::zeros);

  // gradient and shuffled index
  vec grad = vec(tn);
  uvec idx = linspace<uvec>(0, ecn - 1, ecn);

  // active function
  std::shared_ptr<AFmeasure> afgrad = AFfactory().createAFGradient(attrs, arguments);
  std::shared_ptr<AFmeasure> afc = AFfactory().createAFCounts(attrs, arguments);

  for (uword iter = 0; iter < epochs; ++iter) {

    // std::cout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LL(afc->AFCounts(w), efflen, ec, count) << std::endl;
    idx = shuffle(idx);
    uword biter = 0;
    double etai = eta / (1 + decay * iter);

    // mini-batch
    while (biter < ecn) {
      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      // NAG
      grad = afgrad->AFGradient(w - velocity * V, efflen, ec, count, eachidx);
      eg2 = gamma * eg2 + (1 - gamma) * grad % grad;

      // update V
      V = velocity * V + etai / sqrt(eg2 + epsilon) % grad;
      w -= V;

      biter += batchsize;
    }
  }

  // small est & no ec transcripts --> zero
  vec est = afc->AFCounts(w) * cn;
  est.elem(find(est < countLimit)).zeros();
  Rcout << "The log likelihood is " << std::setprecision (20) << LL(est, efflen, ec, count) << "." << std::endl;

  return est;
}


// [[Rcpp::export]]
arma::vec AMSGrad(const arma::vec& efflenraw,
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
  vec m = vec(tn, fill::zeros);
  vec v = vec(tn, fill::zeros);
  uword t = 0;
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
      v = Max(v, beta2 * v + (1 - beta2) * square(grad));
      double etat = eta * sqrt(1 - pow(beta2, t)) / (1 - pow(beta1, t));
      w -= etat * m / (sqrt(v) + epsilon);

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
