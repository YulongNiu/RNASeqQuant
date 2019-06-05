#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::uvec Strsplit(const std::string& s,
                    char delim);

std::vector<arma::uvec> SplitEC(const Rcpp::CharacterVector& ecraw);

std::vector<arma::vec> MatchEfflen(const std::vector<arma::uvec>& ec,
                                   const arma::vec& efflenraw);

arma::vec SpeCount(const arma::vec& est,
                   const arma::uvec& spenumraw);

arma::vec InitAve(const arma::uvec& spenum,
                  const arma::vec& specounts);

arma::vec LambdaSpe(const arma::vec& emlambda,
                    const arma::uvec& spenum,
                    const arma::uvec& spefixcounts);

bool isEqualStr(std::string& str1,
                std::string str2);

#endif
