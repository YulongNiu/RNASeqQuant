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
tmp1 <- EM(plist$efflen, plist$ec, plist$count, 41392)

## ## cpp profiler
## RNASeqEM:::start_profiler("profile.out")
## tmp1 <- RNASeqEM:::EMTest(plist$efflen, plist$ec, plist$count, 41392)
## RNASeqEM:::stop_profiler()

tmp2 <- EMTest(plist$efflen, plist$ec, plist$count, c(20000, 21392))

diffidx <- abs((tmp1 - tmp2))/tmp1 > 0.01
diffidx <- which(diffidx)
cbind(tmp1[diffidx], tmp2[diffidx])

## compare
diffidx <- abs((tmp1 - efflenmat$est_counts))/tmp1 > 0.01
diffidx <- which(diffidx)
cbind(tmp1[diffidx], efflenmat$est_counts[diffidx])




## test share
sourceCpp('../src/TestParal.cpp')

n <- 10
g <- 100000
ecin <- sample(0:9, g*n, replace = TRUE) %>%
  split(1:g)

ecin %>% unlist %>% table

TestShare(ecin, 10)
