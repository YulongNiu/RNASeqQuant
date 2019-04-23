#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include <RcppArmadillo.h>

arma::vec GradientSM_(const arma::vec& w,
                      const std::vector<arma::vec>& efflen,
                      const std::vector<arma::uvec>& ec,
                      const arma::uvec& count,
                      const arma::uvec& idx);

arma::vec GradientSP_(const arma::vec& w,
                      const std::vector<arma::vec>& efflen,
                      const std::vector<arma::uvec>& ec,
                      const arma::uvec& count,
                      const arma::uvec& idx);

arma::vec GradientISRU_(const arma::vec& w,
                        const std::vector<arma::vec>& efflen,
                        const std::vector<arma::uvec>& ec,
                        const arma::uvec& count,
                        const arma::uvec& idx,
                        const double alpha);

#endif
