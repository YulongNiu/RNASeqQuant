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
           const double beta2)
  : tn(tn), beta1(beta1), beta2(beta2){

  t = 0;
  m = vec(tn, fill::zeros);
  v = vec(tn, fill::zeros);

}

// Adam::Adam (const arma::uword tn,
//             const double beta1)
//   : tn(tn), beta1(beta1) {

//   m = vec(tn, fill::zeros);
//   t = 0;
// }


// [[Rcpp::export]]
void Test() {

  Adam testadam(10, 0.9, 0.99);

  Rcout << testadam.tn << endl;
  Rcout << testadam.m << endl;
  Rcout << testadam.v << endl;
  Rcout << testadam.t << endl;
  Rcout << testadam.beta1 << endl;
  Rcout << testadam.beta2 << endl;
}


