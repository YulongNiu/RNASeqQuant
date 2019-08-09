#include "Optfactory.h"
#include "Optimizer.h"
#include "utilities.h"

std::shared_ptr<Optimizer> Optfactory::createOpt(arma::uword tn,
                                                 const Rcpp::List &attrs,
                                                 const Rcpp::List &arguments) {

  std::string optName = attrs["opt"];
  std::shared_ptr<Optimizer> optobj = NULL;

  if (isEqualStr(optName, "Adam")) {
    double beta1 = 0.9;
    double beta2 = 0.999;
    double epsilon = 1e-8;

    if (arguments.containsElementNamed("beta1")) {
      beta1 = arguments["beta1"];
    } else {}

    if (arguments.containsElementNamed("beta2")) {
      beta2 = arguments["beta2"];
    } else {}

    if (arguments.containsElementNamed("epsilon")) {
      beta1 = arguments["epsilon"];
    } else {}

    optobj = std::make_shared<Adam>(tn, beta1, beta2, epsilon);
    // Rcpp::Rcout << *optobj -> beta1 << std::endl;

  }
  else {}

  return optobj;
}

