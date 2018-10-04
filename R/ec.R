##' Read in kallisto output.
##'
##' Read in output files generated from \code{kallisto pseudo} and \code{kallisto quant} (version 0.44.0). The output files include equivalence class (ec), count, and effective length data.
##'
##' @title Standard read in equivalence class (ec) files
##' @param ecpath The path of ec file.
##' @param countpath The path of count file.
##' @param abpath The path of abundance file.
##' @return A \code{list}. 1st element is ec, 2nd element is count of ec, and the 3rd element is effective length of transcript.
##' @examples
##' require('magrittr')
##'
##' ecpath <- system.file('extdata', 'ecexample.ec', package = 'RNASeqEM')
##' countpath <- system.file('extdata', 'count.tsv', package = 'RNASeqEM')
##' abpath <- system.file('extdata', 'abundance.tsv', package = 'RNASeqEM')
##'
##' read_pseudo(ecpath, countpath, abpath)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom utils read.table
##' @importFrom magrittr %>% %<>%
##' @export
##'
read_pseudo <- function(ecpath, countpath, abpath) {

  ec <- read.table(ecpath, stringsAsFactors = FALSE) %>%
    `[`(, 2)

  count <- read.table(countpath, stringsAsFactors = FALSE) %>%
    `[`(, 2)

  efflen <- read.table(abpath, header = TRUE, stringsAsFactors = FALSE) %>%
    `[`(, 3)

  ## remove zero counts
  zeroidx <- count > 0
  count %<>% `[`(zeroidx)
  ec %<>% `[`(zeroidx)

  ecList <- list(ec = ec,
                 count = count,
                 efflen = efflen)

  return(ecList)
}

##' Expectation maximization (EM) model for RNA-seq quantification.
##'
##' EM model for RNA-seq quantification. The equivalence class (ec) with 0 counts are removed, because these counts have no contributes to the final results.
##'
##' @title EM model
##' @param pseudo A \code{list} of ec, count of ec, and effective length of transcripts. It can be the outputs of kallisto.
##' @param maxiter The maximum iteration number with the default value of 10000.
##' @param n The number of CPUs or processors.
##' @inheritParams EMSingle
##' @inheritParams Estcount2Prob
##' @return A \code{numeric vector} indicates probabilities of selecting a read from the different transcripts.
##' @examples
##' ##    f1 f2 f3
##' ## ec1 1 1 1
##' ## ec2 0 1 1
##' ## ec3 1 0 1
##' ## ec4 1 0 0
##' ## ec5 1 1 0
##' plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))
##' EM(plist, spenum = 3, n = 2)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom RcppParallel setThreadOptions
##' @importFrom magrittr %<>%
##' @export
##'
EM <- function(pseudo, spenum, maxiter = 10000, n = 2) {

  ## pseudo information
  count  <- pseudo$count
  ec <- SplitEC(pseudo$ec)
  efflen <-  MatchEfflen(ec, pseudo$efflen)
  spenum %<>% IdxSpenum

  ## set parallel threads
  setThreadOptions(numThreads = n)

  ## stop iteration params from kallisto
  countChangeLimit <- 1e-2
  countChange <- 1e-2
  countLimit <- 1e-8

  ## EM iteration
  startprob <- rep(1/spenum[-1], spenum[-1])
  startcount <- startprob * sum(count) / length(spenum[-1])
  startprob  <- startprob

  for (i in seq_len(maxiter)) {
    est <- EMSingle2(prob = startprob,
                    efflen = efflen,
                    ec = ec,
                    count = count)

    ## stop condition
    stopcond <- est > countChangeLimit &
      (abs(est - startcount)/est) > countChange
    if (sum(stopcond) == 0) {
      ## print(cbind(est, startcount))
      ## print(startprob)
      ## print((abs(est - startcount)/est))
      print(paste0('The iteration number is ', i))
      break
    } else {
      startprob <- Estcount2Prob(est, spenum)
      startcount <- est
    }
  }

  est[est < countLimit] <- 0

  return(est)
}




## ##    f1 f2 f3 f1' f2'
## ## ec1 1  1  0  0  1
## ## ec2 1  0  1  1  0
## ## ec3 0  1  1  0  0
## ## ec4 0  0  0  1  1
## ## ec5 1  0  1  0  1
## ## ec6 1  1  0  0  0

## cp <- c(rep(1/3, 3), rep(1/2, 2))

## cp <- c(c(1/6, 1/2, 1/3), rep(1/2, 2))

## cp <- c(c(1/6, 1/2, 1/3), c(1/4, 3/4))

## effectlen <- list(rep(1, 3), rep(1, 3), rep(1, 2), rep(1, 2), rep(1, 3), rep(1, 2))
## ec <- SplitEC(c('0,1,4', '0,2,3', '1,2', '3,4', '0,2,4', '0,1'))
## effectlen <- MatchEfflen(ec, rep(1, 5))
## ecnum <- rep(1, 6)
## spenum <- IdxSpenum(c(3, 2))

## for (i in 1:100000) {
##   cp <- EMSingle(cp, effectlen, ec, ecnum, spenum)
## }

## cp


library('magrittr')
library('Rcpp')
library('RcppParallel')
library('profvis')
## library('RNASeqEM')
sourceCpp('../src/utilities.cpp')
sourceCpp('../src/EM.cpp')

ecmat <- read.table('/extDisk1/RESEARCH/RNASeqEMtest/athtest/pseudoalignments_ath.tsv', header = TRUE, stringsAsFactors = FALSE)
efflenmat <- read.table('/extDisk1/RESEARCH/RNASeqEMtest/athtest/abundance_ath.tsv', header = TRUE, stringsAsFactors = FALSE)
plist <- list(ec = ecmat$Transcript, count = ecmat$Count, efflen = efflenmat$eff_length)

## cpp
tmp1 <- EMTest(plist$efflen, plist$ec, plist$count, 41392)

## ## cpp profiler
## RNASeqEM:::start_profiler("profile.out")
## tmp1 <- RNASeqEM:::EMTest(plist$efflen, plist$ec, plist$count, 41392)
## RNASeqEM:::stop_profiler()

## R
zeroidx <- ecmat$Count > 0
plist <- list(ec = ecmat$Transcript[zeroidx], count = ecmat$Count[zeroidx], efflen = efflenmat$eff_length)
tmp2 <- EM(plist, 41392)

diffidx <- (tmp1 - tmp2)/tmp1 > 0.01
diffidx <- which(diffidx)
cbind(tmp1[diffidx], tmp2[diffidx])

## compare
diffidx <- (tmp1 - efflenmat$est_counts)/tmp1 > 0.01
diffidx <- which(diffidx)
cbind(tmp1[diffidx], efflenmat$est_counts[diffidx])

## test R version
## profvis(tmp1 <- EM(plist, 41392, n = 8))
## profvis(
##   for (i in 1:2000) {
##     est <- EMSingle(cp, effectlen, ec, ecnum)
##     cp <- Estcount2Prob(est, spenum)
##   }
## )
