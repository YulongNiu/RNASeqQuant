#include <RcppArmadillo.h>

#include "Optimizer.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//======//
// Adam //
//======//
Adam::Adam(const arma::uword tn,
           const double beta1,
           const double beta2,
           double eta,
           const double epsilon)
  : tn(tn), beta1(beta1), beta2(beta2), eta(eta), epsilon(epsilon){

  t = 1;
  m = vec(tn, fill::zeros);
  v = vec(tn, fill::zeros);

}

arma::vec Adam::update (arma::vec& grad,
                        arma::vec& w) {

  m = beta1 * m + (1 - beta1) * grad;
  v = beta2 * v + (1 - beta2) * square(grad);
  double etat = eta * sqrt(1 - pow(beta2, t)) / (1 - pow(beta1, t));

  ++t;

  return w -= etat * m / (sqrt(v) + epsilon);

}

// Adam::Adam (const arma::uword tn,
//             const double beta1)
//   : tn(tn), beta1(beta1) {

//   m = vec(tn, fill::zeros);
//   t = 0;
// }


// [[Rcpp::export]]
void Test() {

  Adam testadam(10, 0.9, 0.99, 0.1, 1e-8);

  Rcout << testadam.tn << endl;
  Rcout << testadam.m << endl;
  Rcout << testadam.v << endl;
  Rcout << testadam.t << endl;
  Rcout << testadam.beta1 << endl;
  Rcout << testadam.beta2 << endl;

  vec g = vec(10); g.fill(0.1);
  vec w = vec(10); w.fill(20);

  Rcout << w << endl;
  Rcout << testadam.update(g, w) << endl;
  Rcout << testadam.t << endl;

}


