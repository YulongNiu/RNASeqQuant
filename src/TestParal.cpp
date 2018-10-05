#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

struct TestShareMem : public Worker
{
  const std::vector<arma::uvec>& ec;
  arma::uvec& estcount;

  TestShareMem(const std::vector<arma::uvec>& ec,
               arma::uvec& estcount)
    : ec(ec), estcount(estcount) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      estcount.elem(ec[i]) += 1;
    }
  }
};


// [[Rcpp::export]]
arma::uvec TestShare(const std::vector<arma::uvec>& ec,
                     const arma::uword num) {

  uvec estcount(num, fill::zeros);

  // create parallel worker and call
  TestShareMem testShareMem(ec, estcount);
  parallelFor(0, ec.size(), testShareMem);

  return estcount;
}

