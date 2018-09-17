#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

struct SumRoot : public Worker
{
  std::vector<arma::vec> idx;
  arma::vec& output;

  SumRoot(std::vector<arma::vec> idx,
          arma::vec& output)
    : idx(idx), output(output) {}

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t i = begin; i < end; ++i) {
      output[i] = sum(idx[i]);
    }
  }
};


// [[Rcpp::export]]
arma::vec SumCppPara(const Rcpp::List x) {

  // allocate the output vector
  vec output(x.size());

  vector<vec> idx(x.size());
  for (int i = 0; i < x.size(); ++i) {
    NumericVector eachelem = x[i];
    vec eachvec(eachelem.begin(), eachelem.size(), false);
    idx[i] = eachvec;
  }

  // SumRoot functor (pass input and output vector)
  SumRoot sumRoot(idx, output);

  // call parallelFor to do the work
  parallelFor(0, x.size(), sumRoot);

  // return the output matrix
  return output;
}

