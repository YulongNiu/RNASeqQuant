#ifndef _SOFTPLUS_H_
#define _SOFTPLUS_H_

#include <RcppArmadillo.h>

arma::vec Logistic(const arma::vec& x);

arma::vec Softplus1(const arma::vec& x);

arma::vec Softplus(const arma::vec& x,
                   const arma::vec& weight);

arma::vec SoftplusGrad1(const arma::vec& x);

arma::vec SoftplusGrad1(const arma::vec& x,
                        const arma::vec& weight);

#endif
