#ifndef AFMEASURE_H_
#define AFMEASURE_H_

#include <RcppArmadillo.h>

using namespace std;
using namespace arma;

typedef arma::vec (*funcPtr)(const arma::vec& w,
                             const std::vector<arma::vec>& efflen,
                             const std::vector<arma::uvec>& ec,
                             const arma::uvec& count,
                             const arma::uvec& idx);

class AFmeasure {
public:
  virtual ~AFmeasure() {};
  virtual arma::vec AFGradient(const arma::vec& w,
                               const std::vector<arma::vec>& efflen,
                               const std::vector<arma::uvec>& ec,
                               const arma::uvec& count,
                               const arma::uvec& idx) = 0;
};

#endif
