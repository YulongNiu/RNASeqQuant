#ifndef _ISRU_H_
#define _ISRU_H_

#include <RcppArmadillo.h>

arma::vec InvSqrtRoot(const arma::vec& x,
                      const double alpha);

arma::vec ISRU1(const arma::vec& x,
                const arma::vec& isr,
                const double alpha);

arma::vec ISRU(const arma::vec& x,
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

#endif
