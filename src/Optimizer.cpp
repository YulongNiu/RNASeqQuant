#include <RcppArmadillo.h>

#include "Optimizer.h"
#include "Optfactory.h"

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


// TestObj(10, list(opt = 'Adam'), list())
// [[Rcpp::export]]
void TestObj(arma::uword tn,
             const Rcpp::List attrs,
             const Rcpp::List arguments) {

  std::shared_ptr<Optimizer> optobj = Optfactory().createOpt(tn, attrs, arguments);
  std::shared_ptr<Optimizer> optobjtest = std::make_shared<Adam>(10, 0.9, 0.999, 1e-8);
  std::shared_ptr<Optimizer> optobjtest2(new Adam(10, 0.9, 0.999, 1e-8));

  vec g = vec(10); g.fill(0.1);
  vec w = vec(10); w.fill(20);

  Rcout << optobjtest -> update(w, g, 0.1) << endl;
}
