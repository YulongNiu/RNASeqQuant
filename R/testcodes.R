## ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test codes~~~~~~~~~~~~~~~~~~~~~~~
## library('Rcpp')
## library('magrittr')
## sourceCpp('../src/utilities.cpp')
## sourceCpp('../src/logsumexp.cpp')
## sourceCpp('../src/likelihood.cpp')
## sourceCpp('../src/GD.cpp')
## sourceCpp('../src/EM.cpp')
## source('ec.R')

## ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~small example~~~~~~~~~~~~~~~~~~~~~~~
## plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))

## w <- c(1, 1, 1)
## Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count, 0:4)
## LL(Softmax1(w), MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)

## EM(plist$efflen, plist$ec, plist$count, spenum = 3)

## Adam(plist$efflen, plist$ec, plist$count, spenum = 3, 300)
## ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~real example~~~~~~~~~~~~~~~~~~~~~
## ## ecpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/athtest/testpseudo/pseudoalignments.ec'
## ## countpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/athtest/testpseudo/pseudoalignments.tsv'
## ## abpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/athtest/testquant/abundance.tsv'
## ## plist <- read_pseudo(ecpath, countpath, abpath)

## plist <- list()
## plist$efflen <- read.table('/extDisk1/RESEARCH/RNASeqQuantTest/GD/abundance_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 3]
## plist$ec <- read.table('/extDisk1/RESEARCH/RNASeqQuantTest/GD/pseudoalignments_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 3]
## plist$count <- read.table('/extDisk1/RESEARCH/RNASeqQuantTest/GD/pseudoalignments_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 2]
## zeroidx <- plist$count > 0
## plist$count %<>% `[`(zeroidx)
## plist$ec %<>% `[`(zeroidx)

## tmp1 <- EM(plist$efflen, plist$ec, plist$count, 41392)
## tmp2 <- Adam(plist$efflen, plist$ec, plist$count, 41392, 500, 1000, 0.01)

## set.seed(12345)
## w <- rnorm(41392, 0, sqrt(1/41392))
## LL(Softmax1(w), MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)
## w <- w - 0.01 * Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count, 0:30000)

## idx <- 0:99
## ## w <- rep(1, 41392)
## w <- rnorm(41392, 0, sqrt(1/41392))
## w <- w - 0.01 * Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count, idx)[idx]
## ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## RNASeqEM:::start_profiler("profile.out")
## tmp1 <- RNASeqEM:::LogSumExpRatio1(rnorm(1:10000))
## RNASeqEM:::stop_profiler()
