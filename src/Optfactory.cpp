#include "Optfactory.h"
#include "Optimizer.h"
#include "utilities.h"

std::shared_ptr<Optimizer> Optfactory::createOpt(arma::uword tn,
                                                 const Rcpp::List &attrs,
                                                 const Rcpp::List &arguments) {

  std::string optName = attrs["opt"];
  std::shared_ptr<Optimizer> optobj = NULL;

  if (isEqualStr(optName, "Adam")) {

    double beta1 = arguments.containsElementNamed("beta1") ? arguments["beta1"] : 0.9;
    double beta2 = arguments.containsElementNamed("beta2") ? arguments["beta1"] : 0.999;
    double epsilon = arguments.containsElementNamed("epsilon") ? arguments["beta1"] : 1e-8;
    optobj = std::make_shared<Adam>(tn, beta1, beta2, epsilon);

  }
  else if (isEqualStr(optName, "NRMSProp")) {

    double gamma = arguments.containsElementNamed("gamma") ? arguments["beta1"] : 0.9;
    double velocity = arguments.containsElementNamed("velocity") ? arguments["beta1"] : 0.9;
    double epsilon = arguments.containsElementNamed("epsilon") ? arguments["beta1"] : 1e-8;
    optobj = std::make_shared<NRMSProp>(tn, gamma, velocity, epsilon);

  }
  else if (isEqualStr(optName, "Adadelta")) {

    double gamma = arguments.containsElementNamed("gamma") ? arguments["beta1"] : 0.9;
    double epsilon = arguments.containsElementNamed("epsilon") ? arguments["beta1"] : 1e-8;
    optobj = std::make_shared<Adadelta>(tn, gamma, epsilon);

  } else {}

  return optobj;
}

