#ifndef AFGRADIENT_H_
#define AFGRADIENT_H_

#include "AFmeasure.h"
#include "gradient.h"
#include "softmax.h"
#include "softplus.h"
#include "isru.h"


//' @inheritParams GradientSM_
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
    return GradientSM_(w, efflen, ec, count, idx);
  }

  arma::vec AFCounts(const arma::vec& w) {
    return Softmax1(w);
  }
};

//' @inheritParams GradientSM_
//' @rdname gradient
//' @keywords internal
//=========//
// Softplus //
//=========//
class AFSP : public AFmeasure {
public:
  arma::vec AFGradient(const arma::vec& w,
                       const std::vector<arma::vec>& efflen,
                       const std::vector<arma::uvec>& ec,
                       const arma::uvec& count,
                       const arma::uvec& idx) {
    return GradientSP_(w, efflen, ec, count, idx);
  }

  arma::vec AFCounts(const arma::vec& w) {
    return Softplus1(w) / sum(Softplus1(w));
  }
};

//' @inheritParams GradientSM_
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
    return GradientISRU_(w, efflen, ec, count, idx, this->alpha);
  }

  arma::vec AFCounts(const arma::vec& w) {
    return ISRU1(w, InvSqrtRoot(w, this->alpha), this->alpha) / sum(ISRU1(w, InvSqrtRoot(w, this->alpha), this->alpha));
  }
};


// //===================================
// // Custom gradient of active function
// //===================================
// class AFCustom : public AFmeasure {
// private:
//   funcGradientPtr fungrad;
//   funcCountsPtr func;
// public:
//   explicit AFCustom(funcGradientPtr fungrad,
//                     funcCountsPtr func) {
//     this->fungrad = fungrad;
//     this->func = func;
//   }
//   ~AFCustom() {}

//   arma::vec AFGradient(const arma::vec& w,
//                        const std::vector<arma::vec>& efflen,
//                        const std::vector<arma::uvec>& ec,
//                        const arma::uvec& count,
//                        const arma::uvec& idx) {
//     return fungrad(w, efflen, ec, count, idx);
//   }

//   arma::vec AFCounts(const arma::vec& w) {
//     return func(w);
//   }
// };

#endif
