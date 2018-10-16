#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


//' Calculate the log-sum-exp calculator
//'
//' \itemize{
//'   \item \code{LogSumExp()}: Weighted log-sum-exp.
//'   \item \code{LogSumExp1()}: log-sum-exp without weight.
//'   \item \code{Softmax()}: Numerator is the exponent of every element of input \code{x}, and denominator is the sum of \code{exp(x)}.
//' }
//'   \item \code{Softmax1()}: \code{weight} is 1.
//' }
//'
//' @title LogSumExp
//' @return
//' \itemize{
//'   \item \code{LogSumExp()} and \code{LogSumExp1(): A \code{double} indicating log-sum-exp.
//'   \item \code{Softmax()}: A \code{arma::vec} number indicate the exp(x_i * weight_i) / sum(exp(x_i * weight_i)).
//'   \item \code{Softmax1()}: A \code{arma::vec} number indicate the exp(x_i) / sum(exp(x_i)).
//' }
//' @param x A \code{arma::vec}.
//' @param weight A \code{arma::vec} indicating the weight.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname logsumexp
//' @keywords internal
// [[Rcpp::export]]
double LogSumExp(const arma::vec& x,
                 const arma::vec& weight) {

  double maxx = max(x);
  vec res = exp(x - maxx) % weight;

  return maxx + log(sum(res));
}


//' @inheritParams LogSumExp
//' @rdname logsumexp
//' @keywords internal
// [[Rcpp::export]]
double LogSumExp1(const arma::vec& x) {

  double maxx = max(x);
  vec res = exp(x - maxx);

  return maxx + log(sum(res));
}


//' @inheritParams LogSumExp
//' @rdname logsumexp
//' @keywords internal
// [[Rcpp::export]]
arma::vec Softmax(const arma::vec& x,
                  const arma::vec& weight) {

  return exp(log(weight) + x - LogSumExp(x, weight));

}


//' @inheritParams LogSumExp
//' @rdname logsumexp
//' @keywords internal
// [[Rcpp::export]]
arma::vec Softmax1(const arma::vec& x) {

  return exp(x - LogSumExp1(x));

}

