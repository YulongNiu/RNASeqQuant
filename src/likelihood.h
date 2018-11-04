#ifndef _LIKELIHOOD_H_
#define _LIKELIHOOD_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

double LLEM(const arma::vec& prob,
            const std::vector<arma::vec>& efflen,
            const std::vector<arma::uvec>& ec,
            const arma::uvec& count);

double LLGD(const arma::vec& prob,
            const std::vector< std::vector< arma::vec > >& efflen,
            const std::vector< std::vector< arma::uvec > >& ec,
            const arma::uvec& count);

#endif

