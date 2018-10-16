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



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test codes~~~~~~~~~~~~~~~~~~~~~~~
library('Rcpp')
library('magrittr')
sourceCpp('../src/utilities.cpp')
sourceCpp('../src/logsumexp.cpp')
sourceCpp('../src/likelihood.cpp')
sourceCpp('../src/GD.cpp')
sourceCpp('../src/EM.cpp')

plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))

Gradient(rep(1, 3), MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)
BGD(plist$efflen, plist$ec, plist$count, spenum = 3, alpha = 0.5)

EM(plist$efflen, plist$ec, plist$count, spenum = 3)



ecpath <- '/extDisk1/RESEARCH/RNASeqEMtest/athtest/testpseudo/pseudoalignments.ec'
countpath <- '/extDisk1/RESEARCH/RNASeqEMtest/athtest/testpseudo/pseudoalignments.tsv'
abpath  <- '/extDisk1/RESEARCH/RNASeqEMtest/athtest/testquant/abundance.tsv'
plist <- read_pseudo(ecpath, countpath, abpath)
tmp1 <- EM(plist$efflen, plist$ec, plist$count, 41392)

w <- rep(1, 41392)
for (i in 1:1) {
  w <- w - 0.1 * Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)
}

tmp1 <- BGD(plist$efflen, plist$ec, plist$count, spenum = 41392, alpha = 0.5, 4)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
