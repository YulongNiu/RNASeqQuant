#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::uvec Strsplit(const std::string& s,
                    char delim);

std::vector<arma::uvec> SplitEC(const Rcpp::CharacterVector& ecraw);

std::vector<arma::vec> MatchEfflen(const std::vector<arma::uvec>& ec,
                                   const arma::vec& efflenraw);

arma::uvec IdxSpenum(const arma::uvec& spenumraw);

void EC2SpeSg(std::vector< arma::uvec >& ecsg,
              std::vector< arma::vec >& efflensg,
              const std::string& ecsgraw,
              const arma::vec& efflenraw,
              const arma::uvec& spenum);

void EC2Spe(std::vector< std::vector< arma::uvec > >& ec,
            std::vector< std::vector< arma::vec > >& efflen,
            const Rcpp::CharacterVector& ecraw,
            const arma::vec& efflenraw,
            const arma::uvec& spenum);

arma::uvec CmpUvec(const std::vector< arma::uvec >& x);

arma::vec CmpVec(const std::vector< arma::vec >& x);

#endif
