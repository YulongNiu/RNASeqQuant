#include "Optfactory.h"
#include "Optimizer.h"
#include "utilities.h"

std::shared_ptr<Optimizer> Optfactory::createOpt(arma::uword tn,
                                                 const Rcpp::List &attrs,
                                                 const Rcpp::List &arguments) {

  std::string optName = attrs["opt"];
  std::shared_ptr<Optimizer> optobj = NULL;

  double epsilon = arguments.containsElementNamed("epsilon") ? arguments["epsilon"] : 1e-8;

  if (isEqualStr(optName, "Momentum")) {
    double gamma = arguments.containsElementNamed("gamma") ? arguments["gamma"] : 0.9;
    optobj = std::make_shared<Momentum>(tn, gamma, epsilon);
  }
  else if (isEqualStr(optName, "NAG")) {
    double gamma = arguments.containsElementNamed("gamma") ? arguments["gamma"] : 0.9;
    double velocity = arguments.containsElementNamed("velocity") ? arguments["velocity"] : 0.9;
    optobj = std::make_shared<NAG>(tn, gamma, velocity, epsilon);
  }
  else if (isEqualStr(optName, "Adagrad")) {
    optobj = std::make_shared<Adagrad>(tn, epsilon);
  }
  else if (isEqualStr(optName, "NAdagrad")) {
    double velocity = arguments.containsElementNamed("velocity") ? arguments["velocity"] : 0.9;
    optobj = std::make_shared<NAdagrad>(tn, velocity, epsilon);
  }
  else if (isEqualStr(optName, "Adadelta")) {
    double gamma = arguments.containsElementNamed("gamma") ? arguments["gamma"] : 0.9;
    optobj = std::make_shared<Adadelta>(tn, gamma, epsilon);
  }
  else if (isEqualStr(optName, "RMSProp")) {
    double gamma = arguments.containsElementNamed("gamma") ? arguments["gamma"] : 0.9;
    optobj = std::make_shared<RMSProp>(tn, gamma, epsilon);
  }
  else if (isEqualStr(optName, "NRMSProp")) {
    double gamma = arguments.containsElementNamed("gamma") ? arguments["gamma"] : 0.9;
    double velocity = arguments.containsElementNamed("velocity") ? arguments["velocity"] : 0.9;
    optobj = std::make_shared<NRMSProp>(tn, gamma, velocity, epsilon);
  }
  else if (isEqualStr(optName, "Adam")) {
    double beta1 = arguments.containsElementNamed("beta1") ? arguments["beta1"] : 0.9;
    double beta2 = arguments.containsElementNamed("beta2") ? arguments["beta2"] : 0.999;
    optobj = std::make_shared<Adam>(tn, beta1, beta2, epsilon);
  }
  else if (isEqualStr(optName, "NAdam")) {
    double beta1 = arguments.containsElementNamed("beta1") ? arguments["beta1"] : 0.9;
    double beta2 = arguments.containsElementNamed("beta2") ? arguments["beta2"] : 0.999;
    optobj = std::make_shared<Adam>(tn, beta1, beta2, epsilon);
  }
  else {}

  return optobj;
}

