#ifndef _ISRU_H_
#define _ISRU_H_

#include <RcppArmadillo.h>

arma::vec InvSqrtRoot(const arma::vec& x,
                      const double alpha);

double ISRU1(const arma::vec& x,
             const arma::vec& isr,
             const double alpha);

double ISRU(const arma::vec& x,
            const arma::vec& isr,
            const arma::vec& weight,
            const double alpha);

arma::vec ISRUGrad1(const arma::vec& x,
                    const arma::vec& isr,
                    const double alpha);

arma::vec ISRUGrad(const arma::vec& x,
                   const arma::vec& isr,
                   const arma::vec& weight,
                   const double alpha);

arma::vec SingleSpeGradISRU(const arma::vec& wlse,
                            const std::vector< arma::vec >& efflensg,
                            const std::vector< arma::vec >& ecd,
                            const arma::vec& ecsum,
                            const arma::vec& ecratio,
                            const arma::uword idx);

arma::vec ECGradISRU(const std::vector< arma::vec >& w,
                     const arma::vec& wlse,
                     const std::vector< arma::vec >& efflensg,
                     const std::vector< arma::vec >& wsg,
                     const double alpha);

#endif
