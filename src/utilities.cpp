#include <RcppArmadillo.h>

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
//' @title Preprocess equivalence classes
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


//' @param ecraw A \code{character vector} and each element is a string with comma delimiter.
//' @rdname strsplit
//' @keywords internal
// [[Rcpp::export]]
std::vector<arma::uvec> SplitEC(const Rcpp::CharacterVector& ecraw) {

  uword ecsize = ecraw.size();
  vector<uvec> res(ecsize);

  for (uword i = 0; i < ecsize; ++i) {
    res[i] = Strsplit(string(ecraw(i)), ',');
  }

  return res;
}


//' Match transcript effect length with equivalence classes.
//'
//' The length of \code{efflen} is equal to number of equivalence classes.
//'
//' @title Match transcript effect length
//' @return A \code{std::vector<arma::vec>} with the same length of \code{ecvec}.
//' @param ec A \code{std::vector<arma::uvec>} containing separated vectors, such as the output of \code{SplitEC()} in this package.
//' @param efflenraw A code{arma::vec} indicating the effect length of transcript.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
std::vector<arma::vec> MatchEfflen(const std::vector<arma::uvec>& ec,
                                   const arma::vec& efflenraw) {
  uword ecsize = ec.size();
  vector<vec> res(ecsize);

  for (uword i = 0; i < ecsize; ++i) {
    res[i] = efflenraw.elem(ec[i]);
  }

  return res;
}


//' Estimated counts of input species.
//'
//' The indices of \code{est} should be consistent with the \code{spenumraw}. For example, the \code{est} is \code{1.5, 2, 3} and \code{spenumraw} is \code{2, 1}, so \code{2.5, 3} will be returned.
//'
//' @title Species estimated counts
//' @return A \code{arma::vec} represents the total estimated counts of each species.
//' @param est A \code{arma::vec} estimated counts of each transcripts.
//' @param spenumraw A \code{arma::uvec} indicated the transcript number in each species.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
arma::vec SpeCount(const arma::vec& est,
                   const arma::uvec& spenumraw) {

  uword sn = spenumraw.n_elem;
  vec res(sn, fill::zeros);

  uword start = 0;

  for (uword i = 0; i < sn; ++i) {
    uword end = start + spenumraw(i) - 1;
    res(i) = sum(est.subvec(start, end));
    start = end + 1;
  }

  return(res);
}
