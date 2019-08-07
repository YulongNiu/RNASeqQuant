#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Optimizer {
public:
  virtual ~Optimizer() {};

  virtual arma::vec update(const arma::vec& w,
                           const arma::vec& grad,
                           const double eta) = 0;
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

    m = vec(tn, fill::zeros);
    v = vec(tn, fill::zeros);
    t = 1;

  }

  arma::vec update(const arma::vec& w,
                   const arma::vec& grad,
                   const double eta) {

    m = beta1 * m + (1 - beta1) * grad;
    v = beta2 * v + (1 - beta2) * square(grad);
    double etat = eta * sqrt(1 - pow(beta2, t)) / (1 - pow(beta1, t));
    vec nextw = w - etat * m / (sqrt(v) + epsilon);
    ++t;

    return nextw;

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

    eg2 = vec(tn, fill::zeros);
    v = vec(tn, fill::zeros);

  }

  // w: estimation
  // grad: gradient
  // eta: learning rate in each epoch
  arma::vec update(const arma::vec& w,
                   const arma::vec& grad,
                   const double eta) {

    eg2 = gamma * eg2 + (1 - gamma) * grad % grad;
    v = velocity * v + eta / sqrt(eg2 + epsilon) % grad;
    vec nextw = w - v;

    return nextw;

  }
};

#endif
