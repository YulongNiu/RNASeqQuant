// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Estcount2Prob
arma::vec Estcount2Prob(const arma::vec& estcount, const arma::uvec& spenum);
RcppExport SEXP _RNASeqQuant_Estcount2Prob(SEXP estcountSEXP, SEXP spenumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type estcount(estcountSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    rcpp_result_gen = Rcpp::wrap(Estcount2Prob(estcount, spenum));
    return rcpp_result_gen;
END_RCPP
}
// EM
arma::vec EM(const arma::vec& efflenraw, const Rcpp::CharacterVector& ecraw, const arma::uvec& countraw, const arma::uvec& spenumraw, const arma::uword maxiter, const arma::uword miniter);
RcppExport SEXP _RNASeqQuant_EM(SEXP efflenrawSEXP, SEXP ecrawSEXP, SEXP countrawSEXP, SEXP spenumrawSEXP, SEXP maxiterSEXP, SEXP miniterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ecraw(ecrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type countraw(countrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenumraw(spenumrawSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type miniter(miniterSEXP);
    rcpp_result_gen = Rcpp::wrap(EM(efflenraw, ecraw, countraw, spenumraw, maxiter, miniter));
    return rcpp_result_gen;
END_RCPP
}
// Adam
arma::vec Adam(const arma::vec& efflenraw, const Rcpp::CharacterVector& ecraw, const arma::uvec& countraw, const arma::uvec& spenumraw, const arma::uword epochs, const arma::uword batchsize, const double alpha);
RcppExport SEXP _RNASeqQuant_Adam(SEXP efflenrawSEXP, SEXP ecrawSEXP, SEXP countrawSEXP, SEXP spenumrawSEXP, SEXP epochsSEXP, SEXP batchsizeSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ecraw(ecrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type countraw(countrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenumraw(spenumrawSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type batchsize(batchsizeSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(Adam(efflenraw, ecraw, countraw, spenumraw, epochs, batchsize, alpha));
    return rcpp_result_gen;
END_RCPP
}
// GradientSM
arma::vec GradientSM(const arma::vec& w, const std::vector<arma::vec>& efflen, const std::vector<arma::uvec>& ec, const arma::uvec& count, const arma::uvec& idx);
RcppExport SEXP _RNASeqQuant_GradientSM(SEXP wSEXP, SEXP efflenSEXP, SEXP ecSEXP, SEXP countSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type efflen(efflenSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::uvec>& >::type ec(ecSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type count(countSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(GradientSM(w, efflen, ec, count, idx));
    return rcpp_result_gen;
END_RCPP
}
// GradientSM2
arma::vec GradientSM2(const arma::vec& w, const std::vector< std::vector< arma::vec > >& efflen, const std::vector< std::vector< arma::uvec > >& ec, const arma::uvec& count, const arma::uvec& spenum, const arma::uvec& idx);
RcppExport SEXP _RNASeqQuant_GradientSM2(SEXP wSEXP, SEXP efflenSEXP, SEXP ecSEXP, SEXP countSEXP, SEXP spenumSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector< arma::vec > >& >::type efflen(efflenSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector< arma::uvec > >& >::type ec(ecSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type count(countSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(GradientSM2(w, efflen, ec, count, spenum, idx));
    return rcpp_result_gen;
END_RCPP
}
// GradientSP
arma::vec GradientSP(const arma::vec& w, const std::vector<arma::vec>& efflen, const std::vector<arma::uvec>& ec, const arma::uvec& count, const arma::uvec& idx);
RcppExport SEXP _RNASeqQuant_GradientSP(SEXP wSEXP, SEXP efflenSEXP, SEXP ecSEXP, SEXP countSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type efflen(efflenSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::uvec>& >::type ec(ecSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type count(countSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(GradientSP(w, efflen, ec, count, idx));
    return rcpp_result_gen;
END_RCPP
}
// LLEM
double LLEM(const arma::vec& prob, const std::vector<arma::vec>& efflen, const std::vector<arma::uvec>& ec, const arma::uvec& count);
RcppExport SEXP _RNASeqQuant_LLEM(SEXP probSEXP, SEXP efflenSEXP, SEXP ecSEXP, SEXP countSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type efflen(efflenSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::uvec>& >::type ec(ecSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type count(countSEXP);
    rcpp_result_gen = Rcpp::wrap(LLEM(prob, efflen, ec, count));
    return rcpp_result_gen;
END_RCPP
}
// LLGD
double LLGD(const arma::vec& prob, const std::vector< std::vector< arma::vec > >& efflen, const std::vector< std::vector< arma::uvec > >& ec, const arma::uvec& count);
RcppExport SEXP _RNASeqQuant_LLGD(SEXP probSEXP, SEXP efflenSEXP, SEXP ecSEXP, SEXP countSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector< arma::vec > >& >::type efflen(efflenSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector< arma::uvec > >& >::type ec(ecSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type count(countSEXP);
    rcpp_result_gen = Rcpp::wrap(LLGD(prob, efflen, ec, count));
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
// SingleSpeGradSM
arma::vec SingleSpeGradSM(const std::vector<arma::vec>& w, const std::vector< arma::vec >& efflensg, const std::vector< arma::vec >& wsg, const arma::vec& ecratio, const arma::uword idx);
RcppExport SEXP _RNASeqQuant_SingleSpeGradSM(SEXP wSEXP, SEXP efflensgSEXP, SEXP wsgSEXP, SEXP ecratioSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<arma::vec>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const std::vector< arma::vec >& >::type efflensg(efflensgSEXP);
    Rcpp::traits::input_parameter< const std::vector< arma::vec >& >::type wsg(wsgSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ecratio(ecratioSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(SingleSpeGradSM(w, efflensg, wsg, ecratio, idx));
    return rcpp_result_gen;
END_RCPP
}
// ECGradSM
arma::vec ECGradSM(const std::vector< arma::vec >& w, const arma::vec& wlse, const std::vector< arma::vec >& efflensg, const std::vector< arma::vec >& wsg);
RcppExport SEXP _RNASeqQuant_ECGradSM(SEXP wSEXP, SEXP wlseSEXP, SEXP efflensgSEXP, SEXP wsgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< arma::vec >& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type wlse(wlseSEXP);
    Rcpp::traits::input_parameter< const std::vector< arma::vec >& >::type efflensg(efflensgSEXP);
    Rcpp::traits::input_parameter< const std::vector< arma::vec >& >::type wsg(wsgSEXP);
    rcpp_result_gen = Rcpp::wrap(ECGradSM(w, wlse, efflensg, wsg));
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
// Test
arma::vec Test(const Rcpp::CharacterVector& ecraw, const arma::vec& efflenraw, const arma::uvec& spenum, const arma::vec& w, const arma::uvec& count, const arma::uvec& idx);
RcppExport SEXP _RNASeqQuant_Test(SEXP ecrawSEXP, SEXP efflenrawSEXP, SEXP spenumSEXP, SEXP wSEXP, SEXP countSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ecraw(ecrawSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type count(countSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(Test(ecraw, efflenraw, spenum, w, count, idx));
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
// IdxSpenum
arma::uvec IdxSpenum(const arma::uvec& spenumraw);
RcppExport SEXP _RNASeqQuant_IdxSpenum(SEXP spenumrawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenumraw(spenumrawSEXP);
    rcpp_result_gen = Rcpp::wrap(IdxSpenum(spenumraw));
    return rcpp_result_gen;
END_RCPP
}
// EC2SpeSg
void EC2SpeSg(std::vector< arma::uvec >& ecsg, std::vector< arma::vec >& efflensg, const std::string& ecsgraw, const arma::vec& efflenraw, const arma::uvec& spenum);
RcppExport SEXP _RNASeqQuant_EC2SpeSg(SEXP ecsgSEXP, SEXP efflensgSEXP, SEXP ecsgrawSEXP, SEXP efflenrawSEXP, SEXP spenumSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector< arma::uvec >& >::type ecsg(ecsgSEXP);
    Rcpp::traits::input_parameter< std::vector< arma::vec >& >::type efflensg(efflensgSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type ecsgraw(ecsgrawSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    EC2SpeSg(ecsg, efflensg, ecsgraw, efflenraw, spenum);
    return R_NilValue;
END_RCPP
}
// EC2Spe
void EC2Spe(std::vector< std::vector< arma::uvec > >& ec, std::vector< std::vector< arma::vec > >& efflen, const Rcpp::CharacterVector& ecraw, const arma::vec& efflenraw, const arma::uvec& spenum);
RcppExport SEXP _RNASeqQuant_EC2Spe(SEXP ecSEXP, SEXP efflenSEXP, SEXP ecrawSEXP, SEXP efflenrawSEXP, SEXP spenumSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::uvec > >& >::type ec(ecSEXP);
    Rcpp::traits::input_parameter< std::vector< std::vector< arma::vec > >& >::type efflen(efflenSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ecraw(ecrawSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type efflenraw(efflenrawSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type spenum(spenumSEXP);
    EC2Spe(ec, efflen, ecraw, efflenraw, spenum);
    return R_NilValue;
END_RCPP
}
// CmpUvec
arma::uvec CmpUvec(const std::vector< arma::uvec >& x);
RcppExport SEXP _RNASeqQuant_CmpUvec(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< arma::uvec >& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(CmpUvec(x));
    return rcpp_result_gen;
END_RCPP
}
// CmpVec
arma::vec CmpVec(const std::vector< arma::vec >& x);
RcppExport SEXP _RNASeqQuant_CmpVec(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector< arma::vec >& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(CmpVec(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RNASeqQuant_Estcount2Prob", (DL_FUNC) &_RNASeqQuant_Estcount2Prob, 2},
    {"_RNASeqQuant_EM", (DL_FUNC) &_RNASeqQuant_EM, 6},
    {"_RNASeqQuant_Adam", (DL_FUNC) &_RNASeqQuant_Adam, 7},
    {"_RNASeqQuant_GradientSM", (DL_FUNC) &_RNASeqQuant_GradientSM, 5},
    {"_RNASeqQuant_GradientSM2", (DL_FUNC) &_RNASeqQuant_GradientSM2, 6},
    {"_RNASeqQuant_GradientSP", (DL_FUNC) &_RNASeqQuant_GradientSP, 5},
    {"_RNASeqQuant_LLEM", (DL_FUNC) &_RNASeqQuant_LLEM, 4},
    {"_RNASeqQuant_LLGD", (DL_FUNC) &_RNASeqQuant_LLGD, 4},
    {"_RNASeqQuant_start_profiler", (DL_FUNC) &_RNASeqQuant_start_profiler, 1},
    {"_RNASeqQuant_stop_profiler", (DL_FUNC) &_RNASeqQuant_stop_profiler, 0},
    {"_RNASeqQuant_LogSumExp", (DL_FUNC) &_RNASeqQuant_LogSumExp, 2},
    {"_RNASeqQuant_LogSumExp1", (DL_FUNC) &_RNASeqQuant_LogSumExp1, 1},
    {"_RNASeqQuant_Softmax", (DL_FUNC) &_RNASeqQuant_Softmax, 2},
    {"_RNASeqQuant_Softmax1", (DL_FUNC) &_RNASeqQuant_Softmax1, 1},
    {"_RNASeqQuant_SingleSpeGradSM", (DL_FUNC) &_RNASeqQuant_SingleSpeGradSM, 5},
    {"_RNASeqQuant_ECGradSM", (DL_FUNC) &_RNASeqQuant_ECGradSM, 4},
    {"_RNASeqQuant_Logistic", (DL_FUNC) &_RNASeqQuant_Logistic, 1},
    {"_RNASeqQuant_Softplus1", (DL_FUNC) &_RNASeqQuant_Softplus1, 1},
    {"_RNASeqQuant_Softplus", (DL_FUNC) &_RNASeqQuant_Softplus, 2},
    {"_RNASeqQuant_SoftplusGrad1", (DL_FUNC) &_RNASeqQuant_SoftplusGrad1, 1},
    {"_RNASeqQuant_SoftplusGrad", (DL_FUNC) &_RNASeqQuant_SoftplusGrad, 2},
    {"_RNASeqQuant_Test", (DL_FUNC) &_RNASeqQuant_Test, 6},
    {"_RNASeqQuant_Strsplit", (DL_FUNC) &_RNASeqQuant_Strsplit, 2},
    {"_RNASeqQuant_SplitEC", (DL_FUNC) &_RNASeqQuant_SplitEC, 1},
    {"_RNASeqQuant_MatchEfflen", (DL_FUNC) &_RNASeqQuant_MatchEfflen, 2},
    {"_RNASeqQuant_IdxSpenum", (DL_FUNC) &_RNASeqQuant_IdxSpenum, 1},
    {"_RNASeqQuant_EC2SpeSg", (DL_FUNC) &_RNASeqQuant_EC2SpeSg, 5},
    {"_RNASeqQuant_EC2Spe", (DL_FUNC) &_RNASeqQuant_EC2Spe, 5},
    {"_RNASeqQuant_CmpUvec", (DL_FUNC) &_RNASeqQuant_CmpUvec, 1},
    {"_RNASeqQuant_CmpVec", (DL_FUNC) &_RNASeqQuant_CmpVec, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_RNASeqQuant(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
