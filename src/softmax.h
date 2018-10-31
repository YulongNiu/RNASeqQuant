#ifndef _SOFTMAX_H_
#define _SOFTMAX_H_

#include <RcppArmadillo.h>

double LogSumExp(const arma::vec& x,
                 const arma::vec& weight);

double LogSumExp1(const arma::vec& x);

arma::vec Softmax(const arma::vec& x,
                  const arma::vec& weight);

arma::vec Softmax1(const arma::vec& x);

arma::vec SingleSpeGradSM(const std::vector<arma::vec>& wnewEfflen,
                          const std::vector<arma::vec>& wnew,
                          const std::vector<arma::vec>& ecEfflen,
                          const std::vector<arma::vec>& ecw,
                          const std::vector<double>& wratio,
                          const arma::uword idx);

arma::vec ECGradSM(const std::vector<arma::vec>& w,
                   const arma::vec wlse,
                   const arma::vec& efflensg,
                   const arma::uvec& ecsg,
                   const arma::uvec& spenum);

#endif
