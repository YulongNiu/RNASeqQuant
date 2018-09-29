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
  const std::vector<arma::vec>& efflen;
  const std::vector<arma::uvec>& es;
  arma::uvec& ecnum;
  arma::vec& newrho;

  ExpectEC(arma::vec& rho,
           const std::vector<arma::vec>& efflen,
           const std::vector<arma::uvec>& es,
           arma::uvec& ecnum,
           arma::vec& newrho)
    : rho(rho), efflen(efflen), es(es), ecnum(ecnum), newrho(newrho) {}

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t i = begin; i < end; ++i) {
      vec eachrho = rho.elem(es[i]) * ecnum[i] / efflen[i];
      newrho.elem(es[i]) += eachrho / sum(eachrho);
    }

  }
};


// [[Rcpp::export]]
arma::vec EMSingle(arma::vec& rho,
                   const std::vector<arma::vec>& efflen,
                   const std::vector<arma::uvec>& es,
                   arma::uvec& ecnum,
                   arma::uvec& spenum) {


  vec newrho(rho.n_elem, fill::zeros);

  // step1: create parallel worker and call
  ExpectEC expectEC(rho, efflen, es, ecnum, newrho);
  parallelFor(0, efflen.size(), expectEC);

  // step2: calculate real newrho
  for (uword i = 0; i < spenum.n_elem - 1; ++i) {
    uword start = spenum(i);
    uword end = spenum(i) + spenum(i+1) - 1;
    newrho.subvec(start, end) /= sum(newrho.subvec(start, end));
  }

  return newrho;
}

