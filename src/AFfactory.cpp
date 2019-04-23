#include "AFfactory.h"
#include "AFgradient.h"
#include "utilities.h"

std::shared_ptr<AFmeasure> AFfactory::createAFGradient(const Rcpp::List &attrs,
                                                       const Rcpp::List &arguments) {

  std::string afName = attrs["method"];
  std::shared_ptr<AFmeasure> afgrad = NULL;

  if (isEqualStr(afName, "Softmax")) {
    afgrad = std::make_shared<GradientSM>();
  }
  else if (isEqualStr(afName, "Softplus")) {
    afgrad = std::make_shared<GradientSP>();
  }
  else if (isEqualStr(afName, "ISRU")) {
    uword alpha = 0.01;
    if (arguments.containsElementNamed("alpha")) {
      alpha = arguments["alpha"];
    } else {}
    afgrad = std::make_shared<GradientISRU>(alpha);
  }
  else if (isEqualStr(afName, "AFCustom")) {
    SEXP func_=arguments["func"];
    funcPtr func = *Rcpp::XPtr<funcPtr>(func_);
    afgrad = std::make_shared<AFCustom>(func);
  }
  else {}

  return afgrad;
}
