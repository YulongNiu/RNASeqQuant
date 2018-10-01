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
  arma::vec& countper;
  const std::vector<arma::vec>& efflen;
  const std::vector<arma::uvec>& ec;
  arma::uvec& ecnum;
  arma::vec& newcountper;

  ExpectEC(arma::vec& countper,
           const std::vector<arma::vec>& efflen,
           const std::vector<arma::uvec>& ec,
           arma::uvec& ecnum,
           arma::vec& newcountper)
    : countper(countper), efflen(efflen), ec(ec), ecnum(ecnum), newcountper(newcountper) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      vec eachcp = countper.elem(ec[i]) / efflen[i];
      newcountper.elem(ec[i]) += eachcp * ecnum(i)/ sum(eachcp);
    }

  }
};


//' Parallel a single EM iteration.
//'
//' Add zero at the head of input \code{spenum}.
//'
//' @title Single EM iteration
//' @return A updated \code{arma::vec} of count percentage.
//' @param countper A \code{arma::vec} of input count percentage.
//' @param efflen A \code{std::vector<arma::vec>} indicated the effective length of transcripts.
//' @param ec A \code{std::vector<arma::vec>} indicated equivalence classes with the same length of \code{efflen}.
//' @param ecnum A \code{arma::uvec} indicated the mapped count number of equivalence class with the same length of \code{efflen}.
//' @param spenum A \code{arma::uvec} indicated number of transcripts.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
arma::vec EMSingle(arma::vec& countper,
                   const std::vector<arma::vec>& efflen,
                   const std::vector<arma::uvec>& ec,
                   arma::uvec& ecnum,
                   arma::uvec& spenum) {

  vec newcountper(countper.n_elem, fill::zeros);

  // step1: create parallel worker and call
  ExpectEC expectEC(countper, efflen, ec, ecnum, newcountper);
  parallelFor(0, efflen.size(), expectEC);

  // step2: calculate real newcountper
  for (uword i = 0; i < spenum.n_elem - 1; ++i) {
    uword start = spenum(i);
    uword end = spenum(i) + spenum(i+1) - 1;
    newcountper.subvec(start, end) /= sum(newcountper.subvec(start, end));
  }

  return newcountper;
}

