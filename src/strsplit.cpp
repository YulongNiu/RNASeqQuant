#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>

#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <iostream>

#include "strsplit.h"

using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// ref: https://ysonggit.github.io/coding/2014/12/16/split-a-string-using-c.html
// [[Rcpp::export]]
arma::uvec Strsplit(const std::string& s,
                    char delim) {

  stringstream ss(s);
  string item;
  uvec tokens(s.size());
  uword i = 0;

  while (getline(ss, item, delim)) {
    tokens(i) = stoul(item);
    ++i;
  }

  uvec res = tokens.subvec(0, i-1);

  return res;
}


// [[Rcpp::export]]
std::vector<arma::uvec> SplitEC(const Rcpp::CharacterVector& ec) {

  uword ecsize = ec.size();
  vector<uvec> res(ecsize);

  for (uword i = 0; i < ecsize; ++i) {
    res[i] = Strsplit(string(ec(i)), ',');
  }

  return res;
}

