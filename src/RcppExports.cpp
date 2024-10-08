// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// EM
Rcpp::List EM(const arma::vec& efflenraw, const Rcpp::CharacterVector& ecraw, const arma::uvec& countraw, const arma::uvec& spenum, const arma::uword maxiter, const arma::uword miniter, const bool details);
RcppExport SEXP _RNASeqQuant_EM(SEXP efflenrawSEXP, SEXP ecrawSEXP, SEXP countrawSEXP, SEXP spenumSEXP, SEXP maxiterSEXP, SEXP miniterSEXP, SEXP detailsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ecraw(ecrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type countraw(countrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type miniter(miniterSEXP);
    Rcpp::traits::input_parameter< const bool >::type details(detailsSEXP);
    rcpp_result_gen = Rcpp::wrap(EM(efflenraw, ecraw, countraw, spenum, maxiter, miniter, details));
    return rcpp_result_gen;
END_RCPP
}
// EMSpe
Rcpp::List EMSpe(const arma::vec& efflenraw, const Rcpp::CharacterVector& ecraw, const arma::uvec& countraw, const arma::uvec& spenum, const arma::vec& spefixcounts, const arma::uword maxiter, const arma::uword miniter, const bool details);
RcppExport SEXP _RNASeqQuant_EMSpe(SEXP efflenrawSEXP, SEXP ecrawSEXP, SEXP countrawSEXP, SEXP spenumSEXP, SEXP spefixcountsSEXP, SEXP maxiterSEXP, SEXP miniterSEXP, SEXP detailsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ecraw(ecrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type countraw(countrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type spefixcounts(spefixcountsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type miniter(miniterSEXP);
    Rcpp::traits::input_parameter< const bool >::type details(detailsSEXP);
    rcpp_result_gen = Rcpp::wrap(EMSpe(efflenraw, ecraw, countraw, spenum, spefixcounts, maxiter, miniter, details));
    return rcpp_result_gen;
END_RCPP
}
// GD
Rcpp::List GD(const arma::vec& efflenraw, const Rcpp::CharacterVector& ecraw, const arma::uvec& countraw, const arma::uvec& spenumraw, const arma::uword epochs, const arma::uword batchsize, const Rcpp::List attrs, const Rcpp::List arguments, const bool details);
RcppExport SEXP _RNASeqQuant_GD(SEXP efflenrawSEXP, SEXP ecrawSEXP, SEXP countrawSEXP, SEXP spenumrawSEXP, SEXP epochsSEXP, SEXP batchsizeSEXP, SEXP attrsSEXP, SEXP argumentsSEXP, SEXP detailsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ecraw(ecrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type countraw(countrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenumraw(spenumrawSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type batchsize(batchsizeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type attrs(attrsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type arguments(argumentsSEXP);
    Rcpp::traits::input_parameter< const bool >::type details(detailsSEXP);
    rcpp_result_gen = Rcpp::wrap(GD(efflenraw, ecraw, countraw, spenumraw, epochs, batchsize, attrs, arguments, details));
    return rcpp_result_gen;
END_RCPP
}
// CountRepeat
arma::uvec CountRepeat(arma::uvec x);
RcppExport SEXP _RNASeqQuant_CountRepeat(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(CountRepeat(x));
    return rcpp_result_gen;
END_RCPP
}
// CountEC
arma::uvec CountEC(const std::vector<arma::uvec>& ec);
RcppExport SEXP _RNASeqQuant_CountEC(SEXP ecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<arma::uvec>& >::type ec(ecSEXP);
    rcpp_result_gen = Rcpp::wrap(CountEC(ec));
    return rcpp_result_gen;
END_RCPP
}
// AdamW
arma::vec AdamW(const arma::vec& efflenraw, const Rcpp::CharacterVector& ecraw, const arma::uvec& countraw, const arma::uvec& ecw, const arma::uvec& spenumraw, const arma::uword epochs, const arma::uword batchsize, const double eta, const Rcpp::List attrs, const Rcpp::List arguments);
RcppExport SEXP _RNASeqQuant_AdamW(SEXP efflenrawSEXP, SEXP ecrawSEXP, SEXP countrawSEXP, SEXP ecwSEXP, SEXP spenumrawSEXP, SEXP epochsSEXP, SEXP batchsizeSEXP, SEXP etaSEXP, SEXP attrsSEXP, SEXP argumentsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ecraw(ecrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type countraw(countrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ecw(ecwSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenumraw(spenumrawSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type batchsize(batchsizeSEXP);
    Rcpp::traits::input_parameter< const double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type attrs(attrsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type arguments(argumentsSEXP);
    rcpp_result_gen = Rcpp::wrap(AdamW(efflenraw, ecraw, countraw, ecw, spenumraw, epochs, batchsize, eta, attrs, arguments));
    return rcpp_result_gen;
END_RCPP
}
// NRMSPropW
arma::vec NRMSPropW(const arma::vec& efflenraw, const Rcpp::CharacterVector& ecraw, const arma::uvec& countraw, const arma::vec& ecw, const arma::uvec& spenumraw, const arma::uword epochs, const arma::uword batchsize, const double eta, const Rcpp::List attrs, const Rcpp::List arguments);
RcppExport SEXP _RNASeqQuant_NRMSPropW(SEXP efflenrawSEXP, SEXP ecrawSEXP, SEXP countrawSEXP, SEXP ecwSEXP, SEXP spenumrawSEXP, SEXP epochsSEXP, SEXP batchsizeSEXP, SEXP etaSEXP, SEXP attrsSEXP, SEXP argumentsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ecraw(ecrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type countraw(countrawSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ecw(ecwSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenumraw(spenumrawSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type batchsize(batchsizeSEXP);
    Rcpp::traits::input_parameter< const double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type attrs(attrsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type arguments(argumentsSEXP);
    rcpp_result_gen = Rcpp::wrap(NRMSPropW(efflenraw, ecraw, countraw, ecw, spenumraw, epochs, batchsize, eta, attrs, arguments));
    return rcpp_result_gen;
END_RCPP
}
// InvSqrtRoot
arma::vec InvSqrtRoot(const arma::vec& x, const double alpha);
RcppExport SEXP _RNASeqQuant_InvSqrtRoot(SEXP xSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(InvSqrtRoot(x, alpha));
    return rcpp_result_gen;
END_RCPP
}
// ISRU1
arma::vec ISRU1(const arma::vec& x, const arma::vec& isr, const double alpha);
RcppExport SEXP _RNASeqQuant_ISRU1(SEXP xSEXP, SEXP isrSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type isr(isrSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(ISRU1(x, isr, alpha));
    return rcpp_result_gen;
END_RCPP
}
// ISRU
arma::vec ISRU(const arma::vec& x, const arma::vec& isr, const arma::vec& weight, const double alpha);
RcppExport SEXP _RNASeqQuant_ISRU(SEXP xSEXP, SEXP isrSEXP, SEXP weightSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type isr(isrSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(ISRU(x, isr, weight, alpha));
    return rcpp_result_gen;
END_RCPP
}
// ISRUGrad1
arma::vec ISRUGrad1(const arma::vec& x, const arma::vec& isr, const double alpha);
RcppExport SEXP _RNASeqQuant_ISRUGrad1(SEXP xSEXP, SEXP isrSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type isr(isrSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(ISRUGrad1(x, isr, alpha));
    return rcpp_result_gen;
END_RCPP
}
// ISRUGrad
arma::vec ISRUGrad(const arma::vec& x, const arma::vec& isr, const arma::vec& weight, const double alpha);
RcppExport SEXP _RNASeqQuant_ISRUGrad(SEXP xSEXP, SEXP isrSEXP, SEXP weightSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type isr(isrSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(ISRUGrad(x, isr, weight, alpha));
    return rcpp_result_gen;
END_RCPP
}
// LL
double LL(const arma::vec& est, const std::vector<arma::vec>& efflen, const std::vector<arma::uvec>& ec, const arma::uvec& count);
RcppExport SEXP _RNASeqQuant_LL(SEXP estSEXP, SEXP efflenSEXP, SEXP ecSEXP, SEXP countSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type est(estSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type efflen(efflenSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::uvec>& >::type ec(ecSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type count(countSEXP);
    rcpp_result_gen = Rcpp::wrap(LL(est, efflen, ec, count));
    return rcpp_result_gen;
END_RCPP
}
// start_profiler
SEXP start_profiler(SEXP str);
RcppExport SEXP _RNASeqQuant_start_profiler(SEXP strSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type str(strSEXP);
    rcpp_result_gen = Rcpp::wrap(start_profiler(str));
    return rcpp_result_gen;
END_RCPP
}
// stop_profiler
SEXP stop_profiler();
RcppExport SEXP _RNASeqQuant_stop_profiler() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(stop_profiler());
    return rcpp_result_gen;
END_RCPP
}
// LogSumExp
double LogSumExp(const arma::vec& x, const arma::vec& weight);
RcppExport SEXP _RNASeqQuant_LogSumExp(SEXP xSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(LogSumExp(x, weight));
    return rcpp_result_gen;
END_RCPP
}
// LogSumExp1
double LogSumExp1(const arma::vec& x);
RcppExport SEXP _RNASeqQuant_LogSumExp1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(LogSumExp1(x));
    return rcpp_result_gen;
END_RCPP
}
// Softmax
arma::vec Softmax(const arma::vec& x, const arma::vec& weight);
RcppExport SEXP _RNASeqQuant_Softmax(SEXP xSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(Softmax(x, weight));
    return rcpp_result_gen;
END_RCPP
}
// Softmax1
arma::vec Softmax1(const arma::vec& x);
RcppExport SEXP _RNASeqQuant_Softmax1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Softmax1(x));
    return rcpp_result_gen;
END_RCPP
}
// Logistic
arma::vec Logistic(const arma::vec& x);
RcppExport SEXP _RNASeqQuant_Logistic(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Logistic(x));
    return rcpp_result_gen;
END_RCPP
}
// Softplus1
arma::vec Softplus1(const arma::vec& x);
RcppExport SEXP _RNASeqQuant_Softplus1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Softplus1(x));
    return rcpp_result_gen;
END_RCPP
}
// Softplus
arma::vec Softplus(const arma::vec& x, const arma::vec& weight);
RcppExport SEXP _RNASeqQuant_Softplus(SEXP xSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(Softplus(x, weight));
    return rcpp_result_gen;
END_RCPP
}
// SoftplusGrad1
arma::vec SoftplusGrad1(const arma::vec& x);
RcppExport SEXP _RNASeqQuant_SoftplusGrad1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(SoftplusGrad1(x));
    return rcpp_result_gen;
END_RCPP
}
// SoftplusGrad
arma::vec SoftplusGrad(const arma::vec& x, const arma::vec& weight);
RcppExport SEXP _RNASeqQuant_SoftplusGrad(SEXP xSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(SoftplusGrad(x, weight));
    return rcpp_result_gen;
END_RCPP
}
// Strsplit
arma::uvec Strsplit(const std::string& s, char delim);
RcppExport SEXP _RNASeqQuant_Strsplit(SEXP sSEXP, SEXP delimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type s(sSEXP);
    Rcpp::traits::input_parameter< char >::type delim(delimSEXP);
    rcpp_result_gen = Rcpp::wrap(Strsplit(s, delim));
    return rcpp_result_gen;
END_RCPP
}
// SplitEC
std::vector<arma::uvec> SplitEC(const Rcpp::CharacterVector& ecraw);
RcppExport SEXP _RNASeqQuant_SplitEC(SEXP ecrawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ecraw(ecrawSEXP);
    rcpp_result_gen = Rcpp::wrap(SplitEC(ecraw));
    return rcpp_result_gen;
END_RCPP
}
// MatchEfflen
std::vector<arma::vec> MatchEfflen(const std::vector<arma::uvec>& ec, const arma::vec& efflenraw);
RcppExport SEXP _RNASeqQuant_MatchEfflen(SEXP ecSEXP, SEXP efflenrawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<arma::uvec>& >::type ec(ecSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    rcpp_result_gen = Rcpp::wrap(MatchEfflen(ec, efflenraw));
    return rcpp_result_gen;
END_RCPP
}
// SpeCount
arma::vec SpeCount(const arma::vec& est, const arma::uvec& spenum);
RcppExport SEXP _RNASeqQuant_SpeCount(SEXP estSEXP, SEXP spenumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type est(estSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    rcpp_result_gen = Rcpp::wrap(SpeCount(est, spenum));
    return rcpp_result_gen;
END_RCPP
}
// InitAve
arma::vec InitAve(const arma::uvec& spenum, const arma::vec& spefixcounts);
RcppExport SEXP _RNASeqQuant_InitAve(SEXP spenumSEXP, SEXP spefixcountsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type spefixcounts(spefixcountsSEXP);
    rcpp_result_gen = Rcpp::wrap(InitAve(spenum, spefixcounts));
    return rcpp_result_gen;
END_RCPP
}
// LambdaSpe
arma::vec LambdaSpe(const arma::vec& emlambda, const arma::uvec& spenum, const arma::vec& spefixcounts);
RcppExport SEXP _RNASeqQuant_LambdaSpe(SEXP emlambdaSEXP, SEXP spenumSEXP, SEXP spefixcountsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type emlambda(emlambdaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type spefixcounts(spefixcountsSEXP);
    rcpp_result_gen = Rcpp::wrap(LambdaSpe(emlambda, spenum, spefixcounts));
    return rcpp_result_gen;
END_RCPP
}
// isEqualStr
bool isEqualStr(std::string& str1, std::string str2);
RcppExport SEXP _RNASeqQuant_isEqualStr(SEXP str1SEXP, SEXP str2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type str1(str1SEXP);
    Rcpp::traits::input_parameter< std::string >::type str2(str2SEXP);
    rcpp_result_gen = Rcpp::wrap(isEqualStr(str1, str2));
    return rcpp_result_gen;
END_RCPP
}
// Max
arma::vec Max(const arma::vec& vec1, const arma::vec& vec2);
RcppExport SEXP _RNASeqQuant_Max(SEXP vec1SEXP, SEXP vec2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type vec2(vec2SEXP);
    rcpp_result_gen = Rcpp::wrap(Max(vec1, vec2));
    return rcpp_result_gen;
END_RCPP
}
// TrueTIdx
arma::uvec TrueTIdx(const std::vector<arma::uvec>& ec);
RcppExport SEXP _RNASeqQuant_TrueTIdx(SEXP ecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<arma::uvec>& >::type ec(ecSEXP);
    rcpp_result_gen = Rcpp::wrap(TrueTIdx(ec));
    return rcpp_result_gen;
END_RCPP
}
// FalseTIdx
arma::uvec FalseTIdx(const std::vector<arma::uvec>& ec, const arma::uvec& spenum);
RcppExport SEXP _RNASeqQuant_FalseTIdx(SEXP ecSEXP, SEXP spenumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<arma::uvec>& >::type ec(ecSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    rcpp_result_gen = Rcpp::wrap(FalseTIdx(ec, spenum));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RNASeqQuant_EM", (DL_FUNC) &_RNASeqQuant_EM, 7},
    {"_RNASeqQuant_EMSpe", (DL_FUNC) &_RNASeqQuant_EMSpe, 8},
    {"_RNASeqQuant_GD", (DL_FUNC) &_RNASeqQuant_GD, 9},
    {"_RNASeqQuant_CountRepeat", (DL_FUNC) &_RNASeqQuant_CountRepeat, 1},
    {"_RNASeqQuant_CountEC", (DL_FUNC) &_RNASeqQuant_CountEC, 1},
    {"_RNASeqQuant_AdamW", (DL_FUNC) &_RNASeqQuant_AdamW, 10},
    {"_RNASeqQuant_NRMSPropW", (DL_FUNC) &_RNASeqQuant_NRMSPropW, 10},
    {"_RNASeqQuant_InvSqrtRoot", (DL_FUNC) &_RNASeqQuant_InvSqrtRoot, 2},
    {"_RNASeqQuant_ISRU1", (DL_FUNC) &_RNASeqQuant_ISRU1, 3},
    {"_RNASeqQuant_ISRU", (DL_FUNC) &_RNASeqQuant_ISRU, 4},
    {"_RNASeqQuant_ISRUGrad1", (DL_FUNC) &_RNASeqQuant_ISRUGrad1, 3},
    {"_RNASeqQuant_ISRUGrad", (DL_FUNC) &_RNASeqQuant_ISRUGrad, 4},
    {"_RNASeqQuant_LL", (DL_FUNC) &_RNASeqQuant_LL, 4},
    {"_RNASeqQuant_start_profiler", (DL_FUNC) &_RNASeqQuant_start_profiler, 1},
    {"_RNASeqQuant_stop_profiler", (DL_FUNC) &_RNASeqQuant_stop_profiler, 0},
    {"_RNASeqQuant_LogSumExp", (DL_FUNC) &_RNASeqQuant_LogSumExp, 2},
    {"_RNASeqQuant_LogSumExp1", (DL_FUNC) &_RNASeqQuant_LogSumExp1, 1},
    {"_RNASeqQuant_Softmax", (DL_FUNC) &_RNASeqQuant_Softmax, 2},
    {"_RNASeqQuant_Softmax1", (DL_FUNC) &_RNASeqQuant_Softmax1, 1},
    {"_RNASeqQuant_Logistic", (DL_FUNC) &_RNASeqQuant_Logistic, 1},
    {"_RNASeqQuant_Softplus1", (DL_FUNC) &_RNASeqQuant_Softplus1, 1},
    {"_RNASeqQuant_Softplus", (DL_FUNC) &_RNASeqQuant_Softplus, 2},
    {"_RNASeqQuant_SoftplusGrad1", (DL_FUNC) &_RNASeqQuant_SoftplusGrad1, 1},
    {"_RNASeqQuant_SoftplusGrad", (DL_FUNC) &_RNASeqQuant_SoftplusGrad, 2},
    {"_RNASeqQuant_Strsplit", (DL_FUNC) &_RNASeqQuant_Strsplit, 2},
    {"_RNASeqQuant_SplitEC", (DL_FUNC) &_RNASeqQuant_SplitEC, 1},
    {"_RNASeqQuant_MatchEfflen", (DL_FUNC) &_RNASeqQuant_MatchEfflen, 2},
    {"_RNASeqQuant_SpeCount", (DL_FUNC) &_RNASeqQuant_SpeCount, 2},
    {"_RNASeqQuant_InitAve", (DL_FUNC) &_RNASeqQuant_InitAve, 2},
    {"_RNASeqQuant_LambdaSpe", (DL_FUNC) &_RNASeqQuant_LambdaSpe, 3},
    {"_RNASeqQuant_isEqualStr", (DL_FUNC) &_RNASeqQuant_isEqualStr, 2},
    {"_RNASeqQuant_Max", (DL_FUNC) &_RNASeqQuant_Max, 2},
    {"_RNASeqQuant_TrueTIdx", (DL_FUNC) &_RNASeqQuant_TrueTIdx, 1},
    {"_RNASeqQuant_FalseTIdx", (DL_FUNC) &_RNASeqQuant_FalseTIdx, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_RNASeqQuant(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
