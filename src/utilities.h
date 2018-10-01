#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::uvec Strsplit(const std::string& s,
                    char delim);

std::vector<arma::uvec> SplitEC(const Rcpp::CharacterVector& ec);

std::vector<arma::vec> MatchEfflen(const std::vector<arma::uvec>& ecvec,
                                   const arma::vec& efflen);

arma::uvec IdxSpenum(const arma::uvec& spenum);

#endif
