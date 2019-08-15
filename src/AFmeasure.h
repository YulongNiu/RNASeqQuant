#ifndef AFMEASURE_H_
#define AFMEASURE_H_

#include <RcppArmadillo.h>

using namespace std;
using namespace arma;

typedef arma::vec (*funcGradientPtr)(const arma::vec& w,
                                     const std::vector<arma::vec>& efflen,
                                     const std::vector<arma::uvec>& ec,
                                     const arma::uvec& count,
                                     const arma::uvec& idx);

typedef arma::vec (*funcCountsPtr)(const arma::vec& w);

class AFmeasure {
public:

  virtual arma::vec AFGradient(const arma::vec& w,
                               const std::vector<arma::vec>& efflen,
                               const std::vector<arma::uvec>& ec,
                               const arma::uvec& count,
                               const arma::uvec& idx) = 0;

  virtual arma::vec AFCounts(const arma::vec& w) = 0;

  virtual ~AFmeasure() {};

};

#endif
