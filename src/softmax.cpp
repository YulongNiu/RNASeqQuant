#include <RcppArmadillo.h>

#include <vector>

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


// emp <- matrix(ncol = 1, nrow = 0)
// SingleSpeGradSM(list(c(1, 1), 1), list(c(1, 1), 1), list(c(1, 1), 1), c(1, 1), 0)
// SingleSpeGradSM(list(c(1, 1), 1), list(1, 1), list(1, 1), c(0.5, 1), 0)
// SingleSpeGradSM(list(c(1, 1), 1), list(1, emp), list(1, emp), c(0.5, 0), 0)
// SingleSpeGradSM(list(c(1, 1), 1), list(c(1, 1), emp), list(c(1, 1), emp), c(1, 0), 0)
// [[Rcpp::export]]
arma::vec SingleSpeGradSM(const std::vector<arma::vec>& w,
                          const std::vector< arma::vec >& efflensg,
                          const std::vector< arma::vec >& wsg,
                          const arma::vec& wratio,
                          const arma::uword idx) {

  // numerator
  vec nr = wsg[idx] + log(1 / efflensg[idx] + sum(wratio) - wratio(idx));

  // denominator
  vec tw = w[idx];
  uword tn = tw.n_elem;
  uword sn = wsg.size();
  vec resw;
  vec reswEfflen;
  for (uword i = 0; i < sn; ++i) {
    if (i != idx && wsg[i].n_elem > 0) {
      resw = join_cols(resw, tw);
      reswEfflen = join_cols(reswEfflen, vec(tn).fill(wratio(i)));
    } else {}
  }

  vec dnw = join_cols(wsg[idx], resw);
  vec dnweight = join_cols(1 / efflensg[idx], reswEfflen);

  // std::cout << nr << std::endl;
  // std::cout << dnw << std::endl;
  // std::cout << dnweight << std::endl;

  double dn = LogSumExp(dnw, dnweight);

  return exp(nr - dn);
}


// emp <- matrix(ncol = 1, nrow = 0)
// ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(c(1, 1), 1), list(c(1, 1), 1))
// ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(1, 1), list(1, 1))
// ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(1, 1), list(1, 1))
// ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(1, emp), list(1, emp))
// ECGradSM(list(c(1, 1), 1), c(log(2)+1, 1), list(c(1, 1), emp), list(c(1, 1), emp))
// [[Rcpp::export]]
arma::vec ECGradSM(const std::vector< arma::vec >& w,
                   const arma::vec wlse,
                   const std::vector< arma::vec >& efflensg,
                   const std::vector< arma::vec >& wsg) {

  // initialization
  uword sn = wsg.size();
  vec wratio(sn, fill::zeros);

  // split each ec
  for (uword i = 0; i < sn; ++i) {
    if (wsg[i].n_elem > 0) {
      wratio(i) = exp(LogSumExp(wsg[i], 1 / efflensg[i]) - wlse(i));
    } else {}
  }

  // calculate each species
  vec res;
  for (uword i = 0; i < sn; ++i) {
    if (wsg[i].n_elem > 0) {
      res = join_cols(res, SingleSpeGradSM(w, efflensg, wsg, wratio, i));
    }
  }

  return res;
}

