#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Optimizer {
public:
  virtual ~Optimizer() {};

  virtual arma::vec update(arma::vec& w,
                           arma::vec& grad,
                           double eta) = 0;
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
  double eta; // learning rate in each epoch
  const double epsilon; // small value, not change

  Adam(const arma::uword tn,
       const double beta1,
       const double beta2,
       const double epsilon);

  arma::vec update(arma::vec& w,
                   arma::vec& grad,
                   double eta);

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
           const double epsilon);

  // w: estimation
  // grad: gradient
  // eta: learning rate in each epoch
  arma::vec update(arma::vec& w,
                   arma::vec& grad,
                   double eta);

};

#endif
