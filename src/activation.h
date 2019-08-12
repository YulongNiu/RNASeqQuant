#ifndef ACTIVATION_H_
#define ACTIVATION_H_

#include "AFmeasure.h"
#include "softmax.h"
#include "softplus.h"
#include "isru.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

//' Gradient.
//'
//' \itemize{
//'   \item \code{AFSM}: Softmax.
//'   \item \code{AFSP}: Softplus.
//'   \item \code{AFISRU}: ISRU.
//' }
//'
//' @title Calculate gradient
//' @return A \code{arma::vec} indicates gradients.
//' @param w A \code{arma::vec} indicates estimated parameters.
//' @param idx A \code{arma::uvec} indicates the index of \code{w} used for gradient descending.
//' @inheritParams LL
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname gradient
//' @keywords internal
//=========//
// Softmax //
//=========//
class AFSM : public AFmeasure {
public:
  arma::vec AFGradient(const arma::vec& w,
                       const std::vector<arma::vec>& efflen,
                       const std::vector<arma::uvec>& ec,
                       const arma::uvec& count,
                       const arma::uvec& idx) {
    vec grad(w.n_elem, fill::zeros);

    for (auto i : idx) {
      grad.elem(ec[i]) += count(i) * Softmax(w.elem(ec[i]), 1/efflen[i]);
    }

    return sum(count.elem(idx)) * Softmax1(w) - grad;
  }

  arma::vec AFCounts(const arma::vec& w) {
    return Softmax1(w);
  }
};

//' @inheritParams AFSM
//' @rdname gradient
//' @keywords internal
//=========//
// Softplus//
//=========//
class AFSP : public AFmeasure {
public:
  arma::vec AFGradient(const arma::vec& w,
                       const std::vector<arma::vec>& efflen,
                       const std::vector<arma::uvec>& ec,
                       const arma::uvec& count,
                       const arma::uvec& idx) {
    vec grad(w.n_elem, fill::zeros);

    for (auto i : idx) {
      grad.elem(ec[i]) += count(i) * SoftplusGrad(w.elem(ec[i]), 1/efflen[i]);
    }

    return sum(count.elem(idx)) * SoftplusGrad1(w) - grad;
  }

  arma::vec AFCounts(const arma::vec& w) {
    return Softplus1(w) / sum(Softplus1(w));
  }
};

//' @inheritParams AFSM
//' @inheritParams InvSqrtRoot
//' @rdname gradient
//' @keywords internal
//=========//
// ISRU    //
//=========//
class AFISRU : public AFmeasure {
private:
  double alpha;
public:
  explicit AFISRU(double alpha) {
    this->alpha = alpha;
  }
  ~AFISRU() {}

  arma::vec AFGradient(const arma::vec& w,
                       const std::vector<arma::vec>& efflen,
                       const std::vector<arma::uvec>& ec,
                       const arma::uvec& count,
                       const arma::uvec& idx) {
    vec grad(w.n_elem, fill::zeros);

    for (auto i : idx) {
      grad.elem(ec[i]) += count(i) * ISRUGrad(w.elem(ec[i]), InvSqrtRoot(w.elem(ec[i]), alpha), 1/efflen[i], alpha);
    }

    return sum(count.elem(idx)) * ISRUGrad1(w, InvSqrtRoot(w, alpha), alpha) - grad;
  }

  arma::vec AFCounts(const arma::vec& w) {
    return ISRU1(w, InvSqrtRoot(w, this->alpha), this->alpha) / sum(ISRU1(w, InvSqrtRoot(w, this->alpha), this->alpha));
  }
};


//===================================//
// Custom gradient of active function//
//===================================//
class AFCustom : public AFmeasure {
private:
  funcGradientPtr funcGrad;
  funcCountsPtr funcC;
public:
  explicit AFCustom(funcGradientPtr funcGrad,
                    funcCountsPtr funcC) {
    this->funcGrad = funcGrad;
    this->funcC = funcC;
  }
  ~AFCustom() {}

  arma::vec AFGradient(const arma::vec& w,
                       const std::vector<arma::vec>& efflen,
                       const std::vector<arma::uvec>& ec,
                       const arma::uvec& count,
                       const arma::uvec& idx) {
    return funcGrad(w, efflen, ec, count, idx);
  }

  arma::vec AFCounts(const arma::vec& w) {
    return funcC(w);
  }
};

#endif
