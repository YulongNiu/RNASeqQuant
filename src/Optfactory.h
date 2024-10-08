#ifndef OPTFACTORY_H_
#define OPTFACTORY_H_

#include <memory>
#include "Optimizer.h"

class Optfactory {
public:

  std::shared_ptr<Optimizer> createOpt(arma::uword tn,
                                       const Rcpp::List& attrs,
                                       const Rcpp::List& arguments);

};

#endif
