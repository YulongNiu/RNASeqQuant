#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Optimizer {
public:
  virtual ~Optimizer() {};

  virtual arma::vec update(arma::vec& grad,
                           arma::vec& w) = 0;
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

  const double beta1; // para1
  const double beta2; // para2
  double eta; // learning rate in each epoch
  const double epsilon; // small value, not change

  Adam(const arma::uword tn,
       const double beta1,
       const double beta2,
       double eta,
       const double epsilon);


  arma::vec update(arma::vec& grad,
                   arma::vec& w);

};

#endif
