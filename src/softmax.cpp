#include <RcppArmadillo.h>

#include <vector>
#include <numeric>

#include "utilities.h"
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
arma::vec SingleSpeGradSM(const std::vector<arma::vec>& wnewEfflen,
                          const std::vector<arma::vec>& wnew,
                          const std::vector<arma::vec>& ecEfflen,
                          const std::vector<arma::vec>& ecw,
                          const std::vector<double>& wratio,
                          const arma::uword idx) {

  // numerator
  vec nr = ecw[idx] + log(1 / ecEfflen[idx] + accumulate(wratio.begin(), wratio.end(), 0.0) - wratio[idx]);

  // denominator
  vec dnw = join_cols(ecw[idx], CbindVector(wnew, idx));
  vec dnweight = join_cols(ecEfflen[idx], CbindVector(wnewEfflen, idx));
  double dn = LogSumExp(dnw, dnweight);


  return exp(nr - dn);
}



// ECGradSM(list(1:3, 4:5), c(LogSumExp1(1:3), LogSumExp1(4:5)), c(1.1, 2.2), c(0, 2), c(0, 3, 2))
// ECGradSM(list(1:3, 4:5), c(LogSumExp1(1:3), LogSumExp1(4:5)), c(1.1, 2.2, 1.3, 3.4), c(0, 2, 3, 4), c(0, 3, 2))
// [[Rcpp::export]]
arma::vec ECGradSM(const std::vector<arma::vec>& w,
                   const arma::vec wlse,
                   const arma::vec& efflensg,
                   const arma::uvec& ecsg,
                   const arma::uvec& spenum) {

  // initialization
  uword sn = spenum.n_elem - 1;
  vector<vec> ecEfflen;
  vector<vec> ecw;
  vector<vec> wnew;
  vector<double> wratio;

  // split each ec
  for (uword i = 0; i < sn; ++i) {

    uword start = spenum(i);
    uword end = spenum(i) + spenum(i+1) - 1;
    uvec eachidx = find(ecsg >= start && ecsg <= end);

    if (eachidx.n_elem > 0) {
      uvec eachec = ecsg.elem(eachidx);
      vec eachefflen = efflensg.elem(eachidx);
      vec eachw = w[i].elem(eachec - start);
      ecEfflen.push_back(eachefflen);
      ecw.push_back(eachw);
      wnew.push_back(w[i]);
      wratio.push_back(exp(LogSumExp(eachw, 1 / eachefflen) - wlse(i)));
    } else {}
  }
  // size equal: ecEfflen ecw wnew wratio wnewEfflen

  // efflen of wnew
  uword realsn = wnew.size();
  vector<vec> wnewEfflen(realsn);
  for (uword i = 0; i < realsn; ++i) {
    wnewEfflen[i] = vec(wnew[i].n_elem).fill(wratio[i]);
  }

  // calculate each species
  vec res(ecsg.n_elem);
  uword start = 0;
  for (uword i = 0; i < realsn; ++i) {

    uword eachlen = ecw[i].n_elem;
    res.subvec(start, eachlen + start - 1) = SingleSpeGradSM(wnewEfflen, wnew, ecEfflen, ecw, wratio, i);
    start += eachlen;

  }

  return res;
}

