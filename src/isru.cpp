#include <RcppArmadillo.h>

#include <vector>

using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


//' Calculate the Inverse square root unit (ISRU)
//'
//' \itemize{
//'   \item \code{InvSqrtRoot()}: Inverse square root.
//'   \item \code{ISRU()} and \code{ISRU1()}: ISRU with or without weight.
//'   \item \code{ISRUGrad()} and \code{ISRUGrad1()}: Internal functions of Partial derivation.
//' }
//'
//' @title ISRU
//' @return
//' \itemize{
//'   \item \code{InvSqrtRoot()}: A \code{arma::vec} indicates the inverse square root \eqn{\frac{1}{sqrt(1+\alpha x^2)}}.
//'   \item \code{ISRU()} and \code{ISRU1()}: A \code{arma::vec} indicates ISRU with (\eqn{(\frac{x_i}{sqrt(1+\alpha x_i^2)} + \frac{1}{\sqrt(\alpha)}) * weight_i} or without weight (\eqn{\frac{x_i}{sqrt(1+\alpha x_i^2)} + \frac{1}{\sqrt(\alpha)}}).
//'   \item \code{ISRUGrad()} and \code{ISRUGrad1()}: A \code{arma::vec} indicates \eqn{(\frac{1}{sqrt(1+\alpha x_i^2)})^3} or \eqn{(\frac{1}{sqrt(1+\alpha x_i^2)})^3 * weight_i}.
//' }
//' @param alpha \code{double}.
//' @inheritParams LogSumExp
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname ISRU
//' @keywords internal
// [[Rcpp::export]]
arma::vec InvSqrtRoot(const arma::vec& x,
                      const double alpha) {

  return 1 / sqrt(1 + alpha * x % x);

}

//' @param isr \code{arma::vec} indicating the inverse square root unit.
//' @inheritParams InvSqrtRoot
//' @rdname ISRU
//' @keywords internal
// [[Rcpp::export]]
arma::vec ISRU1(const arma::vec& x,
                const arma::vec& isr,
                const double alpha) {

  return isr % x + 1 / sqrt(alpha);

}

//' @inheritParams InvSqrtRoot
//' @inheritParams ISRU1
//' @inheritParams LogSumExp
//' @rdname ISRU
//' @keywords internal
// [[Rcpp::export]]
arma::vec ISRU(const arma::vec& x,
               const arma::vec& isr,
               const arma::vec& weight,
               const double alpha) {

  return  isr % x % weight + weight / sqrt(alpha);

}


//' @inheritParams InvSqrtRoot
//' @inheritParams ISRU1
//' @rdname ISRU
//' @keywords internal
// [[Rcpp::export]]
arma::vec ISRUGrad1(const arma::vec& x,
                    const arma::vec& isr,
                    const double alpha) {

  // numerator
  vec nr = pow(isr, 3);

  // denominator
  double dn = sum(isr % x) + x.n_elem / sqrt(alpha);

  return nr / dn;
}


//' @inheritParams InvSqrtRoot
//' @inheritParams ISRU1
//' @inheritParams LogSumExp
//' @rdname ISRU
//' @keywords internal
// [[Rcpp::export]]
arma::vec ISRUGrad(const arma::vec& x,
                   const arma::vec& isr,
                   const arma::vec& weight,
                   const double alpha) {

  // numerator
  vec nr = pow(isr, 3) % weight;

  // denominator
  double dn = sum(isr % x % weight) + sum(weight) / sqrt(alpha);

  return nr / dn;
}
