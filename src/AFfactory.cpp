#include "AFfactory.h"
#include "activation.h"
#include "utilities.h"

std::shared_ptr<AFmeasure> AFfactory::createAF(const Rcpp::List &attrs,
                                               const Rcpp::List &arguments) {
  std::string afName = attrs["af"];
  std::shared_ptr<AFmeasure> af = NULL;

  if (isEqualStr(afName, "Softmax")) {
    af = std::make_shared<AFSM>();
  }
  else if (isEqualStr(afName, "SoftPlus")) {
    af = std::make_shared<AFSP>();
  }
  else if (isEqualStr(afName, "ISRU")) {
    double alpha = arguments.containsElementNamed("alpha") ? arguments["alpha"] : 0.01;
    af = std::make_shared<AFISRU>(alpha);
  }
  else if (isEqualStr(afName, "custom")) {
    SEXP funcGrad_=arguments["gradientAF"];
    SEXP funcC_=arguments["countsAF"];
    funcGradientPtr funcGrad = *Rcpp::XPtr<funcGradientPtr>(funcGrad_);
    funcCountsPtr funcC = *Rcpp::XPtr<funcCountsPtr>(funcC_);
    af = std::make_shared<AFCustom>(funcGrad, funcC);
  }
  else {}

  return af;
}
