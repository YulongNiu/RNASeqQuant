#include <RcppArmadillo.h>

#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>

#include "utilities.h"

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
//' The indices of \code{est} should be consistent with the \code{spenum}. For example, the \code{est} is \code{1.5, 2, 3} and \code{spenum} is \code{2, 1}, so \code{2.5, 3} will be returned.
//'
//' @title Species estimated counts
//' @return A \code{arma::vec} represents the total estimated counts of each species.
//' @param est A \code{arma::vec} estimated counts of each transcripts.
//' @param spenum A \code{arma::uvec} indicated the transcript number in each species.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
arma::vec SpeCount(const arma::vec& est,
                   const arma::uvec& spenum) {

  uword sn = spenum.n_elem;
  vec res(sn, fill::zeros);

  uword start = 0;

  for (uword i = 0; i < sn; ++i) {
    uword end = start + spenum(i) - 1;
    res(i) = sum(est.subvec(start, end));
    start = end + 1;
  }

  return(res);
}


// [[Rcpp::export]]
arma::vec InitAve(const arma::uvec& spenum,
                  const arma::vec& spefixcounts) {

  //spenum.n_elem == specounts_n_elem is true, equal to the number of species

  uword sn = spenum.n_elem;
  vec initest(sum(spenum));

  uword start = 0;

  for (uword i = 0; i < sn; ++i) {
    uword end = start + spenum(i) - 1;
    initest.subvec(start, end).fill(spefixcounts(i) / spenum(i));
    start = end + 1;
  }

  return(initest);
}


// [[Rcpp::export]]
arma::vec LambdaSpe(const arma::vec& emlambda,
                    const arma::uvec& spenum,
                    const arma::vec& spefixcounts) {

  // emlambda.n_elem is number of transcripts
  // spenum.n_elem == spefixcounts.n_elem is true, species number
  uword sn = spenum.n_elem;
  uword tn = sum(spenum);
  vec res(tn);

  uword start = 0;

  for (uword i = 0; i < sn; ++i) {
    uword end = start + spenum(i) - 1;
    vec eachemlambda = emlambda.subvec(start, end);
    res.subvec(start, end) = eachemlambda * spefixcounts(i) / sum(eachemlambda);
    start = end + 1;
  }

  return(res);
}


//' Compare two strings
//'
//' @title Two strings comparison
//' @return A \code{bool} indicating whether two strings are equal.
//' @param str1 str2 \code{std::string} strings.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
bool isEqualStr(std::string& str1,
                std::string str2) {
  return str1.compare(str2) == 0;
}



//' Pairwise maximum elements of two vectors
//'
//' @title Extract maximum elements
//' @return A \code{arma::vec} vector indicating maximum element by pair-wise comparison.
//' @param vec1 vec2 \code{arma::vec} vectors.
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @keywords internal
// [[Rcpp::export]]
arma::vec Max(const arma::vec& vec1,
              const arma::vec& vec2) {

  mat m = join_rows(vec1, vec2);
  return max(m, 1);

}


// [[Rcpp::export]]
arma::uvec TrueTIdx(const std::vector<arma::uvec>& ec) {

  // collapse ec
  uvec coles;
  for (auto i : ec) {
    coles = join_cols(coles, i);
  }

  return unique(coles);
}


// [[Rcpp::export]]
arma::uvec FalseTIdx(const std::vector<arma::uvec>& ec,
                     const arma::uvec& spenum) {

  uword tn = sum(spenum);
  uvec ftIdx = zeros<uvec>(tn);

  ftIdx.elem(TrueTIdx(ec)).ones();

  return find(ftIdx == 0);
}

// [[Rcpp::export]]
arma::vec PreInit(const std::vector<arma::vec>& efflen,
                  const std::vector<arma::uvec>& ec,
                  const arma::uvec& count,
                  const arma::uvec& spenum) {
  vec init(sum(spenum), fill::zeros);

  for (uword i = 0; i < ec.size(); ++i) {
    vec eachcp = 1 / efflen[i];
    init.elem(ec[i]) += eachcp * count(i) / sum(eachcp);
  }

  vec initnorm = init / sum(init);
  for (uword i = 0; i < initnorm.n_elem; ++i) {
    initnorm[i] = (initnorm[i] == 0) ? -1e8 : log(initnorm[i]);
  }

  return initnorm;
}
