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
           const double epsilon)
  : tn(tn), beta1(beta1), beta2(beta2), epsilon(epsilon){

  m = vec(tn, fill::zeros);
  v = vec(tn, fill::zeros);
  t = 1;

}

arma::vec Adam::update (arma::vec& w,
                        arma::vec& grad,
                        double eta) {

  m = beta1 * m + (1 - beta1) * grad;
  v = beta2 * v + (1 - beta2) * square(grad);
  double etat = eta * sqrt(1 - pow(beta2, t)) / (1 - pow(beta1, t));
  w -= etat * m / (sqrt(v) + epsilon);

  ++t;

  return w;

}


//==========//
// NRMSProp //
//==========//
NRMSProp::NRMSProp(const arma::uword tn,
                   const double gamma,
                   const double velocity,
                   const double epsilon)
  : tn(tn), gamma(gamma), velocity(velocity), epsilon(epsilon){

  eg2 = vec(tn, fill::zeros);
  v = vec(tn, fill::zeros);

}

arma::vec NRMSProp::update (arma::vec& w,
                            arma::vec& grad,
                            double eta) {

  eg2 = gamma * eg2 + (1 - gamma) * grad % grad;
  v = velocity * v + eta / sqrt(eg2 + epsilon) % grad;
  w -= v;

  return w;

}


// [[Rcpp::export]]
void Test() {

  Adam testadam(10, 0.9, 0.99, 1e-8);

  Rcout << testadam.tn << endl;
  Rcout << testadam.m << endl;
  Rcout << testadam.v << endl;
  Rcout << testadam.t << endl;
  Rcout << testadam.beta1 << endl;
  Rcout << testadam.beta2 << endl;

  vec g = vec(10); g.fill(0.1);
  vec w = vec(10); w.fill(20);

  Rcout << w << endl;
  Rcout << testadam.update(w, g, 0.1) << endl;
  Rcout << testadam.t << endl;

}


