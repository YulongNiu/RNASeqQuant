#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

struct ExpectEC : public Worker
{
  arma::vec& rho;
  const std::vector<arma::vec>& effectlen;
  const std::vector<arma::uvec>& ts;
  arma::uvec& ecnum;
  arma::vec& newrho;

  ExpectEC(arma::vec& rho,
           const std::vector<arma::vec>& effectlen,
           const std::vector<arma::uvec>& ts,
           arma::uvec& ecnum,
           arma::vec& newrho)
    : rho(rho), effectlen(effectlen), ts(ts), ecnum(ecnum), newrho(newrho) {}

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t i = begin; i < end; ++i) {
      vec eachrho = rho.elem(ts[i]) * ecnum[i] / effectlen[i];
      newrho.elem(ts[i]) += eachrho / sum(eachrho);
    }

  }
};


// [[Rcpp::export]]
arma::vec EMSingle(arma::vec& rho,
                   const std::vector<arma::vec>& effectlen,
                   const std::vector<arma::uvec>& ts,
                   arma::uvec& ecnum) {

  // allocate the output vector
  vec newrho(rho.n_elem, fill::zeros);

  ExpectEC expectEC(rho, effectlen, ts, ecnum, newrho);

  parallelFor(0, effectlen.size(), expectEC);

  return newrho;
}

