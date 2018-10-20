#ifndef _LOGSUMEXP_H_
#define _LOGSUMEXP_H_

#include <RcppArmadillo.h>

double LogSumExp(const arma::vec& x,
                 const arma::vec& weight);

double LogSumExp1(const arma::vec& x);

arma::vec Softmax(const arma::vec& x,
                  const arma::vec& weight);

arma::vec Softmax1(const arma::vec& x);

#endif
