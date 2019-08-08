#include <RcppArmadillo.h>

#include "Optimizer.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
void Test() {

  Optimizer *testadam  = new NRMSProp(10, 0.9, 0.9, 1e-8);

  vec g = vec(10); g.fill(0.1);
  vec w = vec(10); w.fill(20);

  Rcout << w << endl;
  Rcout << testadam -> update(w, g, 0.1) << endl;
  Rcout << testadam -> preupdate(w) << endl;

}


