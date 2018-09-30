#ifndef _STRSPLIT_H_
#define _STRSPLIT_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::uvec Strsplit(const std::string& s,
                    char delim);

std::vector<arma::uvec> SplitEC(const Rcpp::CharacterVector& ec);

std::vector<arma::vec> MatchEfflen(const std::vector<arma::uvec>& ecvec,
                                   const arma::vec& efflen);

#endif
