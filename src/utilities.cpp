#include <RcppArmadillo.h>
#define ARMA_64BIT_WORD 1

#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


//' Split strings and equivalence classes.
//'
//' \itemize{
//'   \item \code{Strsplit()}: Split a \code{string} with an user-defined delimiter.
//'   \item \code{SplitEC()}: Split batch of equivalence classes in \code{string} format.
//' }
//'
//' @title Preprocess equivalence classes.
//' @return
//' \itemize{
//'   \item \code{Strsplit()}: A \code{arma::uvec} indicating the corresponding transcripts ID (starts from 0).
//'   \item \code{SplitEC}: A \code{std::vector<arma::uvec>} and each element indicates the transcripts IDs.
//' }
//' @param s A \code{string}.
//' @param delim The delimiter
//' @references \href{https://ysonggit.github.io/coding/2014/12/16/split-a-string-using-c.html}{split string in C++}
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @rdname strsplit
//' @keywords internal
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


//' @param ec A \code{character vector} and each element is a string with comma delimiter.
//' @rdname strsplit
//' @keywords internal
// [[Rcpp::export]]
std::vector<arma::uvec> SplitEC(const Rcpp::CharacterVector& ec) {

  uword ecsize = ec.size();
  vector<uvec> res(ecsize);

  for (uword i = 0; i < ecsize; ++i) {
    res[i] = Strsplit(string(ec(i)), ',');
  }

  return res;
}


//' Match transcript effect length with equivalence classes.
//'
//' The length of \code{efflen} is equal to number of transcripts.
//'
//' @title Match transcript effect length
//' @return A \code{std::vector<arma::vec>} with the same length of \code{ecvec}.
//' @param ecvec A \code{std::vector<arma::uvec>} containing separated vectors, such as the output of \code{SplitEC()} in this package.
//' @param efflen A code{arma::vec} indicating the effect length of transcript.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
std::vector<arma::vec> MatchEfflen(const std::vector<arma::uvec>& ecvec,
                                   const arma::vec& efflen) {
  uword ecsize = ecvec.size();
  vector<vec> res(ecsize);

  for (uword i = 0; i < ecsize; ++i) {
    res[i] = efflen.elem(ecvec[i]);
  }

  return res;
}


//' Index transcripts number of input species.
//'
//' Add zero at the head of input \code{spenum}.
//'
//' @title Index species number
//' @return A \code{arma::uvec}.
//' @param spenum A \code{arma::uvec} indicated the transcript number in each species.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
arma::uvec IdxSpenum(const arma::uvec& spenum) {

  uvec res(spenum.n_elem + 1, fill::zeros);

  for (uword i = 0; i < spenum.n_elem; ++i) {
    res(i+1) = spenum(i);
  }

  return res;
}
