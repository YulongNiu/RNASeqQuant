#include "AFfactory.h"
#include "activation.h"
#include "utilities.h"

std::shared_ptr<AFmeasure> AFfactory::createAFGradient(const Rcpp::List &attrs,
                                                       const Rcpp::List &arguments) {

  std::string afName = attrs["af"];
  std::shared_ptr<AFmeasure> afgrad = NULL;

  if (isEqualStr(afName, "Softmax")) {
    afgrad = std::make_shared<AFSM>();
  }
  else if (isEqualStr(afName, "Softplus")) {
    afgrad = std::make_shared<AFSP>();
  }
  else if (isEqualStr(afName, "ISRU")) {
    double alpha = arguments.containsElementNamed("alpha") ? arguments["alpha"] : 0.01;
    afgrad = std::make_shared<AFISRU>(alpha);
  }
  // else if (isEqualStr(afName, "custom")) {
  //   SEXP func_=arguments["gradientAF"];
  //   funcGradientPtr func = *Rcpp::XPtr<funcGradientPtr>(func_);
  //   afgrad = std::make_shared<AFCustom>(func);
  // }
  else {}

  return afgrad;
}


std::shared_ptr<AFmeasure> AFfactory::createAFCounts(const Rcpp::List &attrs,
                                                     const Rcpp::List &arguments) {

  std::string afName = attrs["af"];
  std::shared_ptr<AFmeasure> afc = NULL;

  if (isEqualStr(afName, "Softmax")) {
    afc = std::make_shared<AFSM>();
  }
  else if (isEqualStr(afName, "Softplus")) {
    afc = std::make_shared<AFSP>();
  }
  else if (isEqualStr(afName, "ISRU")) {
    double alpha = arguments.containsElementNamed("alpha") ? arguments["alpha"] : 0.01;
    afc = std::make_shared<AFISRU>(alpha);
  }
  // else if (isEqualStr(afName, "custom")) {
  //   SEXP func_=arguments["countsAF"];
  //   funcCountsPtr func = *Rcpp::XPtr<funcCountsPtr>(func_);
  //   afc = std::make_shared<AFCustom>(func);
  // }
  else {}

  return afc;
}
