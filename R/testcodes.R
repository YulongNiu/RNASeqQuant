##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test codes~~~~~~~~~~~~~~~~~~~~~~~
library('Rcpp')
library('magrittr')
sourceCpp('../src/utilities.cpp')
sourceCpp('../src/softmax.cpp')
sourceCpp('../src/softplus.cpp')
sourceCpp('../src/likelihood.cpp')
sourceCpp('../src/GD.cpp')
sourceCpp('../src/EM.cpp')
source('ec.R')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~small example~~~~~~~~~~~~~~~~~~~~~~~
plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))

w <- c(1, 1, 1)
Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count, 0:4)
LL(Softmax1(w), MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)

EM(plist$efflen, plist$ec, plist$count, spenum = 3)

Adam(plist$efflen, plist$ec, plist$count, spenum = 3, 400)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~real example~~~~~~~~~~~~~~~~~~~~~
## ecpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/athtest/testpseudo/pseudoalignments.ec'
## countpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/athtest/testpseudo/pseudoalignments.tsv'
## abpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/athtest/testquant/abundance.tsv'
## plist <- read_pseudo(ecpath, countpath, abpath)

plist <- list()
plist$efflen <- read.table('/home/Yulong/RESEARCH/RNASeqQuantTest/GD/abundance_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 3]
plist$ec <- read.table('/home/Yulong/RESEARCH/RNASeqQuantTest/GD/pseudoalignments_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 3]
plist$count <- read.table('/home/Yulong/RESEARCH/RNASeqQuantTest/GD/pseudoalignments_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 2]

zeroidx <- plist$count > 0
plist$count %<>% `[`(zeroidx)
plist$ec %<>% `[`(zeroidx)


## kallisto EM
kallistoest <- read.table('/home/Yulong/RESEARCH/RNASeqQuantTest/GD/abundance_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 4]

## RNASeqQuant EM
emest <- EM(plist$efflen, plist$ec, plist$count, 41392)

## RNASeqQuant GD
gdest <- Adam(plist$efflen, plist$ec, plist$count, 41392, 200, 1000, 0.01)

## merge res
mergeres <- cbind(kallistoest, emest, gdest)
colnames(mergeres) <- c('kallistoest', 'emest', 'gdest')

## correlation corefficient
cor(mergeres, method = 'pearson')
cor(mergeres, method = 'spearman')

set.seed(12345)
w <- rnorm(41392, 0, sqrt(1/41392))
LL(Softmax1(w), MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)
w <- w - 0.01 * Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count, 0:30000)

idx <- 0:99
## w <- rep(1, 41392)
w <- rnorm(41392, 0, sqrt(1/41392))
w <- w - 0.01 * Gradient(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count, idx)[idx]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


RNASeqEM:::start_profiler("profile.out")
tmp1 <- RNASeqEM:::LogSumExpRatio1(rnorm(1:10000))
RNASeqEM:::stop_profiler()
