#include <RcppArmadillo.h>

#include <algorithm>
#include <vector>

#include "utilities.h"
#include "likelihood.h"
#include "AFfactory.h"
#include "Optfactory.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


//' Gradient decent (GD) model for RNA-seq quantification.
//'
//' GD model for RNA-seq quantification. The equivalence class (ec) with 0 counts are removed, because these counts have no contributes to the final results.
//'
//' @title GD model
//' @param batchsize Mini-batch size, it should be smaller or equal to \code{epochs}.
//' @param attrs Set active functions and optimization algorithms.
//' \itemize{
//'   \item \code{af}: "Softmax", "SoftPlus", "ISRU".
//'   \item \code{gd}: "Adagrad", "NAdagrad", "Adadelta", "RMSProp", "NRMSProp", "Adam", "NAdam", and "AMSGrad".
//' }
//' @param arguments Set advance parameters for \code{attrs}.
//' \itemize{
//'   \item \code{alpha}: used in "ISRU", default 0.01.
//'   \item \code{eta} and \code{decay}: Learning rate and decay rate, default 0.1 and 0.03, respectively. Used in "Adagrad", "NAdagrad", "RMSProp", "NRMSProp", "Adam", "NAdam", and "AMSGrad".
//'   \item \code{gamma}: used in "Adadelta", "RMSProp", and "NRMSProp", default 0.9.
//'   \item \code{velocity}: used in "NAdagrad" and "NRMSProp", default 0.9.
//'   \item \code{beta1} and \code{beta2}: used in "Adam" , "NAdam", "AdaMax", and "AMSGrad", default 0.9 and 0.999, respectively.
//'   \item \code{epsilon}: used in all optimization algorithms, default 1e-08.
//'   \item \code{assign0}: used in "AdaMax" indicating whether assign 0 to transcripts that having no mapped reads. \code{FALSE} is recommended for "AdaMax".
//'   \item \code{gradientAF} and \code{countsAF}: customized active function for gradient and counts.
//' }
//' @inheritParams EM
//' @return A \code{List} indicates estimated counts of transcripts.
//' @examples
//' ## Single species
//' ##    f1 f2 f3
//' ## ec1 1 1 1
//' ## ec2 0 1 1
//' ## ec3 1 0 1
//' ## ec4 1 0 0
//' ## ec5 1 1 0
//' plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'Softmax', opt = 'Adagrad'), list(eta = 0.5, decay = 0.03))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'Softmax', opt = 'NAdagrad'), list(eta = 0.5, decay = 0.03))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'Softmax', opt = 'Adadelta'), list(gamma = 0.8))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'Softmax', opt = 'RMSProp'), list(eta = 0.1, decay = 0.03))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'Softmax', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'Softmax', opt = 'Adam'), list(eta = 0.1, decay = 0.03))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'Softmax', opt = 'NAdam'), list(eta = 0.1, decay = 0.03))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'Softmax', opt = 'AdaMax'), list(eta = 0.1, decay = 0.03, assign0 = FALSE))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'Softmax', opt = 'AMSGrad'), list(eta = 0.1, decay = 0.03))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'SoftPlus', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03))
//' GD(plist$efflen, plist$ec, plist$count, spenum = 3,
//'    list(af = 'ISRU', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03))
//'
//' ## Two species
//' ##    f1 f2 f3 f1' f2'
//' ## ec1 1  1  0  0  1
//' ## ec2 1  0  1  1  0
//' ## ec3 0  1  1  0  0
//' ## ec4 0  0  0  1  1
//' ## ec5 1  0  1  0  1
//' ## ec6 1  1  0  0  0
//' plist <- list(ec = c('0,1,4', '0,2,3', '1,2', '3,4', '0,2,4', '0,1'),
//'               count = rep(1, 6), efflen = rep(1, 5))
//' GD(plist$efflen, plist$ec, plist$count, spenum = c(3, 2),
//'    list(af = 'Softmax', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03), batchsize = 1024)
//' @author Yulong Niu \email{yulong.niu@@hotmail.com}
//' @export
// [[Rcpp::export]]
Rcpp::List GD(const arma::vec& efflenraw,
              const Rcpp::CharacterVector& ecraw,
              const arma::uvec& countraw,
              const arma::uvec& spenum,
              const Rcpp::List attrs,
              const Rcpp::List arguments,
              const arma::uword maxiter = 10000,
              const arma::uword miniter = 50,
              const arma::uword batchsize = 1024,
              const bool details = false) {

  // stop iteration settings from kallisto
  double countChangeLimit = 1e-2;
  double countChange = 1e-2;
  double countLimit = 1e-8;

  // step1: pseudo information remove zero counts
  uvec zeros = find(countraw > 0);
  IntegerVector zerosidx(zeros.begin(), zeros.end());

  uvec count = countraw.elem(zeros);
  vector<uvec> ec = SplitEC(ecraw[zerosidx]);
  vector<vec> efflen = MatchEfflen(ec, efflenraw);

  // step2: initialization
  uword tn = sum(spenum);
  uword cn = sum(count);
  uword ecn = ec.size();
  vec resll(maxiter, fill::zeros);

  //~~~~~~~~~~~~~~~~~old init~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Glorot normal initializer/Xavier normal initializer
  // vec w(tn); w.fill(0.01);
  // vec w = randn<vec>(tn) / sqrt(tn);
  // uvec ftidx = FalseTIdx(ec, spenum);
  // w.elem(ftidx).fill(-1e8);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  vec w = PreInit(efflen, ec, count, spenum);

  // gd settings
  if (arguments.containsElementNamed("assign0") && !arguments["assign0"]) {
    w.elem(find(w < -1e7)).fill(0.01);
  } else {}
  double eta = arguments.containsElementNamed("eta") ? arguments["eta"] : 0.1;
  double decay = arguments.containsElementNamed("decay") ? arguments["decay"] : 0.03;

  // active function
  std::shared_ptr<AFmeasure> af = AFfactory().createAF(attrs, arguments);
  // GD method
  std::shared_ptr<Optimizer> gd = Optfactory().createOpt(tn, attrs, arguments);

  vec startest = af->AFCounts(w) * cn;
  vec est(tn, fill::zeros);

  // step3: gradient decent
  vec grad = vec(tn);
  uvec idx = linspace<uvec>(0, ecn - 1, ecn);
  uword iter;
  for (iter = 0; iter < maxiter; ++iter) {
    Rcout << std::setprecision (10) << min(w) << "|" << max(w) << "|" << LogSumExp1(w) << "|" << LL(af->AFCounts(w) * cn, efflen, ec, count) << std::endl;

    // record running details
    if (details) {
      resll(iter) = LL(startest, efflen, ec, count);
    } else {}

    idx = shuffle(idx);
    uword biter = 0;
    double etai = eta / (1 + decay * iter);

    // mini-batch
    while (biter < ecn) {

      uword endi = biter + batchsize - 1;
      endi = (endi >= ecn) ? (ecn - 1) : endi;
      uvec eachidx = idx.subvec(biter, endi);

      grad = af->AFGradient(gd->preupdate(w), efflen, ec, count, eachidx);
      w = gd->update(w, grad, etai);

      biter += batchsize;
    }

    est = af->AFCounts(w) * cn;

    uword nopassn = sum((est > countChangeLimit) %
                        (abs(est - startest) > countChange));

    // Rcout << nopassn << endl;

    if (nopassn == 0 && iter >= miniter - 1) {
      break;
    } else {
      startest = est;
    }
  }

  // check if maxiter
  if (iter == maxiter) {--iter;} else {}

  // step4: small est & no ec transcripts --> zero
  est.elem(find(est < countLimit)).zeros();

  List res = List::create(_["counts"] = est,
                          _["ll"] = resll.subvec(0, iter));

  Rcout << "The iteration number is " << iter + 1
        << ". The log likelihood is " << std::setprecision (20) << LL(est, efflen, ec, count)
        << "." << std::endl;

  return res;
}
