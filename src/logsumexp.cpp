#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


//' Calculate the log-sum-exp calculator
//'
//' \itemize{
//'   \item \code{LogSumExp()}: Weighted log-sum-exp.
//'   \item \code{LogSumExpRatio()}: Numerator is the exponent of every element of input \code{x}, and denominator is the sum of \code{exp(x)}.
//' }
//'   \item \code{LogSumExpRatio1()}: \code{weight} is 1.
//' }
//'
//' @title LogSumExp
//' @return
//' \itemize{
//'   \item \code{LogSumExp()}: A \code{double} indicating log-sum-exp.
//'   \item \code{LogSumExpRatio()}: A \code{arma::vec} number indicate the exp(x_i * weight_i) / sum(exp(x_i * weight_i)).
//'   \item \code{LogSumExpRatio1()}: A \code{arma::vec} number indicate the exp(x_1) / sum(exp(x_i)).
//' }
//' @param x A \code{arma::vec}.
//' @param weight A \code{arma::vec} indicating the weight.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname logsumexp
//' @keywords internal
// [[Rcpp::export]]
double LogSumExp(const arma::vec& x,
                 const arma::vec& weight) {

  vec res(x.n_elem, fill::zeros);
  double maxx = max(x);

  res = exp(x - maxx) % weight;

  return maxx + log(sum(res));
}


//' @inheritParams LogSumExp
//' @rdname logsumexp
//' @keywords internal
// [[Rcpp::export]]
arma::vec LogSumExpRatio(const arma::vec& x,
                         const arma::vec& weight) {

  vec res(x.n_elem, fill::zeros);

  for (uword i = 0; i < x.n_elem; ++i) {
    res(i) =  weight(i) / exp(LogSumExp(x - x(i), weight));
  }

  return res;
}


//' @inheritParams LogSumExp
//' @rdname logsumexp
//' @keywords internal
// [[Rcpp::export]]
arma::vec LogSumExpRatio1(const arma::vec& x) {

  vec weight(x.n_elem, fill::ones);

  return LogSumExpRatio(x, weight);
}

