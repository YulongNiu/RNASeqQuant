#ifndef _LOGSUMEXP_H_
#define _LOGSUMEXP_H_

#include <RcppArmadillo.h>

double LogSumExp(const arma::vec& x,
                 const arma::vec& weight);

arma::vec LogSumExpRatio(const arma::vec& x,
                         const arma::vec& weight);

arma::vec LogSumExpRatio1(const arma::vec& x);

#endif
