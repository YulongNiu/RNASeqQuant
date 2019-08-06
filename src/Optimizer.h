#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Optimizer {
public:
  virtual ~Optimizer() {};

  // virtual void update(const arma::vec grad) = 0;
};


//======//
// Adam //
//======//
class Adam : public Optimizer {
public:
  const arma::uword tn; // transcription number
  arma::vec m;
  arma::vec v;
  arma::uword t;
  const double beta1;
  const double beta2;

  Adam(const arma::uword tn,
       const double beta1,
       const double beta2);

  // void init(const arma::uword tn);

  void update(const arma::vec grad);

};

// class Adam : public Optimizer {
// public:
//   const arma::uword tn; // transcription number
//   arma::vec& m;
//   arma::uword t;
//   const double beta1;

//   Adam(const arma::uword tn,
//        const double beta1);

// };


#endif
