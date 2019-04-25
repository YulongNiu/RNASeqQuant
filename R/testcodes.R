##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test codes~~~~~~~~~~~~~~~~~~~~~~~
library('Rcpp')
library('magrittr')
sourceCpp('../src/utilities.cpp')
sourceCpp('../src/softmax.cpp')
sourceCpp('../src/softplus.cpp')
sourceCpp('../src/isru.cpp')
sourceCpp('../src/likelihood.cpp')
sourceCpp('../src/GD.cpp')
sourceCpp('../src/EM.cpp')
source('ec.R')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~small example~~~~~~~~~~~~~~~~~~~~~~~
plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))

EM(plist$efflen, plist$ec, plist$count, spenum = 3) %>% .$counts

Adam(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'Softmax'), list())
Adam(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'Softplus'), list())
Adam(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'ISRU'), list(alpha = 0.1))

Adagrad(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'Softmax'), list())
Adagrad(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'Softplus'), list())
Adagrad(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'ISRU'), list(alpha = 0.1))

RMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'Softmax'), list())
RMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'Softplus'), list())
RMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'ISRU'), list(alpha = 0.1))

NRMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'Softmax'), list())
NRMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'Softplus'), list())
NRMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(method = 'ISRU'), list(alpha = 0.1))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~real example~~~~~~~~~~~~~~~~~~~~~
## ecpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/athtest/testpseudo/pseudoalignments.ec'
## countpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/athtest/testpseudo/pseudoalignments.tsv'
## abpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/athtest/testquant/abundance.tsv'
## plist <- read_pseudo(ecpath, countpath, abpath)

## ecpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/simulate/human/humanpseudo/pseudoalignments.ec'
## countpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/simulate/human/humanpseudo/pseudoalignments.tsv'
## abpath <- '/extDisk1/RESEARCH/RNASeqQuantTest/simulate/human/humanquant/abundance.tsv'
## plist <- read_pseudo(ecpath, countpath, abpath)

plist <- list()
plist$efflen <- read.table('/extDisk1/RESEARCH/RNASeqQuantTest/GD/abundance_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 3]
plist$ec <- read.table('/extDisk1/RESEARCH/RNASeqQuantTest/GD/pseudoalignments_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 3]
plist$count <- read.table('/extDisk1/RESEARCH/RNASeqQuantTest/GD/pseudoalignments_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 2]

zeroidx <- plist$count > 0
plist$count %<>% `[`(zeroidx)
plist$ec %<>% `[`(zeroidx)

## kallisto EM
## kallistoest <- read.table('/extDisk1/RESEARCH/RNASeqQuantTest/simulate/human/humanquant/abundance.tsv', stringsAsFactors = FALSE, header = TRUE)[, 4]
kallistoest <- read.table('/extDisk1/RESEARCH/RNASeqQuantTest/GD/abundance_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 4]

## RNASeqQuant EM
emest <- EM(plist$efflen, plist$ec, plist$count, length(plist$efflen), detail = FALSE)

## RNASeqQuant GD
gdest <- Adam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'Softmax'), list())
gdest <- Adam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'Softplus'), list())
gdest <- Adam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'ISRU'), list(alpha = 0.1))


gdest <- Adagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'Softmax'), list())
gdest <- Adagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'Softplus'), list())
gdest <- Adagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'ISRU'), list(alpha = 0.1))

gdest <- RMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'Softmax'), list())
gdest <- RMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'Softplus'), list())
gdest <- RMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'ISRU'), list(alpha = 0.1))

gdest <- NRMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'Softmax'), list())
gdest <- NRMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'Softplus'), list())
gdest <- NRMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(method = 'ISRU'), list(alpha = 0.1))


## nice test
## Adam mini-batch
gdest <- Adam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 300, 1024, 0.1, list(method = 'Softmax'), list())
gdest <- NAdam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 300, 1024, 0.1, list(method = 'Softmax'), list())
gdest <- NAdagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 300, 1024, 0.1, list(method = 'Softmax'), list())


## NRMSProp full batch
emest <- EM(plist$efflen, plist$ec, plist$count, length(plist$efflen), detail = TRUE)
gdest <- NRMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 300, 36580, 0.005, list(method = 'Softmax'), list())
gdest <- NAdagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 300, 36580, 0.1, list(method = 'Softmax'), list())

## merge res
mergeres <- cbind(kallistoest, emest, gdest)
colnames(mergeres) <- c('kallistoest', 'emest', 'gdest')

## correlation coefficient
cor(mergeres, method = 'pearson')
cor(mergeres, method = 'spearman')

## test full batch
set.seed(12345)
w <- rnorm(41392, 0, sqrt(1/41392))

LL(Softmax1(w), MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)
w <- w - 0.01 * GradientSM(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count, length(plist$efflen), 0:30000)

LL(Softplus1(w)/sum(Softplus1(w)), MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count)
w <- w - 0.01 * GradientSP(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count, 0:30000)

## w <- rnorm(41392, 0, sqrt(1/41392))
w <- w - 0.01 * GradientSM(w, MatchEfflen(SplitEC(plist$ec), plist$efflen), SplitEC(plist$ec), plist$count, idx)[idx+1]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## RNASeqEM:::start_profiler("profile.out")
## tmp1 <- RNASeqEM:::LogSumExpRatio1(rnorm(1:10000))
## RNASeqEM:::stop_profiler()
