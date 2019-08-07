#include <RcppArmadillo.h>

#include "Optimizer.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

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


