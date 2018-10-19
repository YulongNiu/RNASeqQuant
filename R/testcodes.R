## ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test codes~~~~~~~~~~~~~~~~~~~~~~~
## library('Rcpp')
## library('magrittr')
## sourceCpp('../src/utilities.cpp')
## sourceCpp('../src/logsumexp.cpp')
## sourceCpp('../src/likelihood.cpp')
## sourceCpp('../src/GD.cpp')
## sourceCpp('../src/EM.cpp')
## source('ec.R')

## plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))

## w <- c(1, 2, 3)
## Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)
## LL(Softmax1(w), MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)

## BGD(plist$efflen, plist$ec, plist$count, spenum = 3, alpha = 0.5)

## EM(plist$efflen, plist$ec, plist$count, spenum = 3)



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
## LL(rep(1/41392, 41392), MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)
## ## LLR(rep(1/41392, 41392), MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)

## set.seed(12345)
## w <- rnorm(41392, 0, sqrt(1/41392))

## w <- rep(1, 41392)
## for (i in 1:2) {
##   w <- w - 0.01 * Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)
## }

## head(Softmax1(w))
## sum(Softmax1(w) * sum(plist$count))
## head(Softmax1(w) * sum(plist$count))


## idx <- 36000:36099
## w <- rep(1, 41392)
## Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen)[idx], SplitEC(plist$ec)[idx], plist$count[idx])[idx]

## for (i in 1:10) {
##   batchSize <- 10000
##   idx <- sample(1 : length(plist$count), batchSize)
##   w <- w - 0.01 * Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen)[idx], SplitEC(plist$ec)[idx], plist$count[idx])
## }

## head(Softmax1(w))
## sum(Softmax1(w) * sum(plist$count))
## head(Softmax1(w) * sum(plist$count))

## ## est <- Estw2Estcount(startw, sum(plist$count))
## ## prob <- LogSumExpRatio1(w)

## tmp1 <- BGD(plist$efflen, plist$ec, plist$count, spenum = 41392, alpha = 0.5, 2)
## ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## RNASeqEM:::start_profiler("profile.out")
## tmp1 <- RNASeqEM:::LogSumExpRatio1(rnorm(1:10000))
## RNASeqEM:::stop_profiler()
