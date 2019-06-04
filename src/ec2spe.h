#ifndef _EC2SPE_H_
#define _EC2SPE_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::uvec AccuSum(const arma::uvec spenum);

arma::uvec EC2SpeEach(const arma::uvec accuIdx,
                      const arma::uvec ec);

std::vector<arma::uvec> EC2Spe(const arma::uvec accuIdx,
                               const std::vector<arma::uvec> ec);

#endif
