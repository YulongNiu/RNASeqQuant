#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>

using namespace std;

// ref: https://ysonggit.github.io/coding/2014/12/16/split-a-string-using-c.html
// [[Rcpp::export]]
std::vector<unsigned long> Strsplit(const std::string& s,
                                    char delim) {

  stringstream ss(s);
  string item;
  vector<unsigned long> tokens;
  while (getline(ss, item, delim)) {
    tokens.push_back(stoul(item));
  }
  return tokens;
}


