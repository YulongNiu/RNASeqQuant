#ifndef AFFACTORY_H_
#define AFFACTORY_H_

#include <memory>
#include "AFmeasure.h"

//==============================
// AF factory
//==============================
class AFfactory {
public:
  std::shared_ptr<AFmeasure> createAF(const Rcpp::List& attrs,
                                      const Rcpp::List& arguments);
};

#endif
