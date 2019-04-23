#ifndef AFGRADIENT_H_
#define AFGRADIENT_H_

#include "AFmeasure.h"
#include "gradient.h"


//' @inheritParams GradientSM_
//' @rdname gradient
//' @keywords internal
//=========//
// Softmax //
//=========//
class GradientSM : public AFmeasure {
public:
  arma::vec AFGradient(const arma::vec& w,
                       const std::vector<arma::vec>& efflen,
                       const std::vector<arma::uvec>& ec,
                       const arma::uvec& count,
                       const arma::uvec& idx) {
    return GradientSM_(w, efflen, ec, count, idx);
  }
};

//' @inheritParams GradientSM_
//' @rdname gradient
//' @keywords internal
//=========//
// Softplus //
//=========//
class GradientSP : public AFmeasure {
public:
  arma::vec AFGradient(const arma::vec& w,
                       const std::vector<arma::vec>& efflen,
                       const std::vector<arma::uvec>& ec,
                       const arma::uvec& count,
                       const arma::uvec& idx) {
    return GradientSP_(w, efflen, ec, count, idx);
  }
};

//' @inheritParams GradientSM_
//' @inheritParams InvSqrtRoot
//' @rdname gradient
//' @keywords internal
//=========//
// ISRU    //
//=========//
class GradientISRU : public AFmeasure {
private:
  double alpha;
public:
  explicit GradientISRU(double alpha) {
    this->alpha = alpha;
  }
  ~GradientISRU() {}
  arma::vec AFGradient(const arma::vec& w,
                       const std::vector<arma::vec>& efflen,
                       const std::vector<arma::uvec>& ec,
                       const arma::uvec& count,
                       const arma::uvec& idx) {
    return GradientISRU_(w, efflen, ec, count, idx, this->alpha);
  }
};


//===================================
// Custom gradient of active function
//===================================
class AFCustom : public AFmeasure {
private:
  funcPtr func;
public:
  explicit AFCustom(funcPtr func) {
    this->func = func;
  }
  ~AFCustom() {}
  arma::vec AFGradient(const arma::vec& w,
                       const std::vector<arma::vec>& efflen,
                       const std::vector<arma::uvec>& ec,
                       const arma::uvec& count,
                       const arma::uvec& idx) {
    return func(w, efflen, ec, count, idx);
  }
};

#endif
