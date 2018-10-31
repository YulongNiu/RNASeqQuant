#include <RcppArmadillo.h>

#include <vector>

#include "softmax.h"

using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


//' Calculate the log-sum-exp and softmax calculator
//'
//' \itemize{
//'   \item \code{LogSumExp()}: Weighted log-sum-exp.
//'   \item \code{LogSumExp1()}: log-sum-exp without weight.
//'   \item \code{Softmax()}: Numerator is the exponent of every element of input \code{x}, and denominator is the sum of \code{exp(x)}.
//'   \item \code{Softmax1()}: \code{weight} is 1.
//' }
//'
//' @title Softmax
//' @return
//' \itemize{
//'   \item \code{LogSumExp()} and \code{LogSumExp1()}: A \code{double} indicates log-sum-exp.
//'   \item \code{Softmax()}: A \code{arma::vec} indicates the \eqn{\frac{\mathrm{e}^{x_i} * weight_i}{\sum{\mathrm{e}^{x_i} * weight_i}}}.
//'   \item \code{Softmax1()}: A \code{arma::vec} indicates the \eqn{\frac{\mathrm{e}^{x_i}}{\sum{\mathrm{e}^{x_i}}}}.
//' }
//' @param x A \code{arma::vec}.
//' @param weight A \code{arma::vec} indicating the weight.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname softmax
//' @keywords internal
// [[Rcpp::export]]
double LogSumExp(const arma::vec& x,
                 const arma::vec& weight) {

  double maxx = max(x);

  return maxx + log(sum(exp(x - maxx) % weight));
}


//' @inheritParams LogSumExp
//' @rdname softmax
//' @keywords internal
// [[Rcpp::export]]
double LogSumExp1(const arma::vec& x) {

  double maxx = max(x);

  return maxx + log(sum(exp(x - maxx)));
}


//' @inheritParams LogSumExp
//' @rdname softmax
//' @keywords internal
// [[Rcpp::export]]
arma::vec Softmax(const arma::vec& x,
                  const arma::vec& weight) {

  return exp(log(weight) + x - LogSumExp(x, weight));

}


//' @inheritParams LogSumExp
//' @rdname softmax
//' @keywords internal
// [[Rcpp::export]]
arma::vec Softmax1(const arma::vec& x) {

  return exp(x - LogSumExp1(x));

}


// [[Rcpp::export]]
arma::vec EachGradSM(const std::vector<arma::vec>& wsplit,
                     const std::vector<arma::vec>& efflen,
                     const std::vector<arma::vec>& w,
                     const arma::vec wratio,
                     const uword idx) {
  // denominator

}


// [[Rcpp::export]]
arma::vec CutSM(const std::vector<arma::vec>& wsplit,
                const arma::vec& efflensg,
                const arma::uvec& ecsg,
                const arma::uvec& spenum) {

  // initialization
  uword sn = spenum.n_elem - 1;
  vector<vec> efflen;
  vector<vec> w;
  vector<vec> wsplitnew;
  vec wratio(sn, fill::zeros);

  // split each ec
  for (uword i = 0; i < sn; ++i) {

    uword start = spenum(i);
    uword end = spenum(i) + spenum(i+1) - 1;
    uvec eachidx = find(ecsg >= start && ecsg <= end);

    if (eachidx.n_elem > 0) {
      uvec eachec = ecsg.elem(eachidx);
      vec eachefflen = efflensg.elem(eachidx);
      vec eachw = wsplit[i].elem(eachec - start);
      efflen.push_back(eachefflen);
      w.push_back(eachw);
      wsplitnew.push_back(wsplit[i]);
      wratio(i) = exp(LogSumExp(eachw, 1 / eachefflen) - LogSumExp1(wsplit[i]));
    } else {}

  }

  for (auto s : w) {
    std::cout << s << std::endl;
  }

  for (auto s : efflen) {
    std::cout << s << std::endl;
  }

  for (auto s : wsplitnew) {
    std::cout << s << std::endl;
  }

  std::cout << wratio.elem(find(wratio != 0)) << std::endl;

  return efflensg;
}

