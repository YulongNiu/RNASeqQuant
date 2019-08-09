#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include <RcppArmadillo.h>
#include <cmath>

class Optimizer {
public:
  virtual ~Optimizer() {};

  // w: estimation
  // grad: gradient
  // eta: learning rate in each epoch
  virtual arma::vec update(const arma::vec& w,
                           const arma::vec& grad,
                           const double eta) = 0;

  virtual arma::vec preupdate(const arma::vec& w) = 0;

};


//======//
// Adam //
//======//
class Adam : public Optimizer {
public:
  const arma::uword tn; // #transcripts

  arma::vec m;
  arma::vec v;
  arma::uword t;

  const double beta1; // para1 for Adam
  const double beta2; // para2 for Adam
  const double epsilon; // small value, not change

  Adam(const arma::uword tn,
       const double beta1,
       const double beta2,
       const double epsilon)
    : tn(tn), beta1(beta1), beta2(beta2), epsilon(epsilon) {

    m = arma::vec(tn, arma::fill::zeros);
    v = arma::vec(tn, arma::fill::zeros);
    t = 1;

  }

  arma::vec update(const arma::vec& w,
                   const arma::vec& grad,
                   const double eta) {

    m = beta1 * m + (1 - beta1) * grad;
    v = beta2 * v + (1 - beta2) * arma::square(grad);
    double etat = eta * std::sqrt(1 - std::pow(beta2, t)) / (1 - std::pow(beta1, t));
    arma::vec nextw = w - etat * m / (arma::sqrt(v) + epsilon);
    ++t;

    return nextw;

  }

  arma::vec preupdate(const arma::vec& w) {
    return w;
  }
};


//==========//
// NRMSProp //
//==========//
class NRMSProp : public Optimizer {
public:
  const arma::uword tn; // #transcripts

  arma::vec eg2;
  arma::vec v;

  const double gamma; // para1 for RMSProp
  const double velocity; // para2 for NAG
  const double epsilon; // small value, not change

  NRMSProp(const arma::uword tn,
           const double gamma,
           const double velocity,
           const double epsilon)
    : tn(tn), gamma(gamma), velocity(velocity), epsilon(epsilon) {

    eg2 = arma::vec(tn, arma::fill::zeros);
    v = arma::vec(tn, arma::fill::zeros);

  }

  arma::vec update(const arma::vec& w,
                   const arma::vec& grad,
                   const double eta) {

    eg2 = gamma * eg2 + (1 - gamma) * grad % grad;
    v = velocity * v + eta / arma::sqrt(eg2 + epsilon) % grad;
    arma::vec nextw = w - v;

    return nextw;

  }

  arma::vec preupdate(const arma::vec& w) {
    return w - velocity * v;
  }
};


//===========//
// Adadelta //
//=========//
class Adadelta : public Optimizer {
public:
  const arma::uword tn; // #transcripts

  arma::vec eg2;
  arma::vec edx2;

  const double gamma; // para1 for Adadelta
  const double epsilon; // small value, not change

  Adadelta(const arma::uword tn,
           const double gamma,
           const double epsilon)
    : tn(tn), gamma(gamma), epsilon(epsilon) {

    eg2 = arma::vec(tn, arma::fill::zeros);
    edx2 = arma::vec(tn, arma::fill::zeros);

  }

  arma::vec update(const arma::vec& w,
                   const arma::vec& grad,
                   const double eta) {

    eg2 = gamma * eg2 + (1 - gamma) * grad % grad;
    arma::vec dx = -arma::sqrt(edx2 + epsilon) / arma::sqrt(eg2 + epsilon) % grad;
    edx2 = gamma * edx2 + (1 - gamma) * dx % dx;
    arma::vec nextw = w + dx;

    return nextw;

  }

  arma::vec preupdate(const arma::vec& w) {
    return w;
  }
};

#endif


