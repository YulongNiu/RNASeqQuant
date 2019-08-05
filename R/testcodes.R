##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test codes~~~~~~~~~~~~~~~~~~~~~~~
library('tibble')
library('dplyr')
library('ggplot2')
library('Rcpp')
library('magrittr')
sourceCpp('../src/utilities.cpp')
sourceCpp('../src/softmax.cpp')
sourceCpp('../src/softplus.cpp')
sourceCpp('../src/isru.cpp')
sourceCpp('../src/likelihood.cpp')
sourceCpp('../src/GD.cpp')
sourceCpp('../src/WGD.cpp')
sourceCpp('../src/EM.cpp')
source('ec.R')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~small example~~~~~~~~~~~~~~~~~~~~~~~
plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))

EM(plist$efflen, plist$ec, plist$count, spenum = 3) %>% .$counts

Adam(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'Softmax'), list())
Adam(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'Softplus'), list())
Adam(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'ISRU'), list(alpha = 0.1))

Adagrad(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'Softmax'), list())
Adagrad(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'Softplus'), list())
Adagrad(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'ISRU'), list(alpha = 0.1))

RMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'Softmax'), list())
RMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'Softplus'), list())
RMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'ISRU'), list(alpha = 0.1))

NRMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'Softmax'), list())
NRMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'Softplus'), list())
NRMSProp(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, 0.1, list(af = 'ISRU'), list(alpha = 0.1))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~real example~~~~~~~~~~~~~~~~~~~~~
testpath <-'/extDisk1/RESEARCH/RNASeqQuantTestPython/athtest'
## testpath <- '/extDisk1/RESEARCH/RNASeqQuantTestPython/simulate/human/humanpseudo'
ecpath <- file.path(testpath, 'testpseudo/pseudoalignments.ec')
countpath <- file.path(testpath, 'testpseudo/pseudoalignments.tsv')
abpath <- file.path(testpath, 'testquant/abundance.tsv')
plist <- read_pseudo(ecpath, countpath, abpath)
w <- read.table(ecpath, stringsAsFactors = FALSE) %>%
  `[`(, 2) %>%
  SplitEC %>%
  CountEC

## plist <- list()
## plist$efflen <- read.table('abundance_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 3]
## plist$ec <- read.table('/extDisk1/RESEARCH/RNASeqQuantTestPython/GD/pseudoalignments_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 3]
## plist$count <- read.table('/extDisk1/RESEARCH/RNASeqQuantTestPython/GD/pseudoalignments_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 2]
## zeroidx <- plist$count > 0
## plist$count %<>% `[`(zeroidx)
## plist$ec %<>% `[`(zeroidx)

## kallisto EM
## kallistoest <- read.table('/extDisk1/RESEARCH/RNASeqQuantTestPython/simulate/human/humanquant/abundance.tsv', stringsAsFactors = FALSE, header = TRUE)[, 4]
kallistoest <- read.table('/extDisk1/RESEARCH/RNASeqQuantTestPython/GD/abundance_ath.tsv', stringsAsFactors = FALSE, header = TRUE)[, 4]

## RNASeqQuant EM
emest <- EMSpe(plist$efflen, plist$ec, plist$count, c(10000, 10000, 21392), rep(3731388, 3), detail = TRUE)
emest <- EM(plist$efflen, plist$ec, plist$count, c(10000, 10000, 21392), detail = TRUE)
emest <- EM(plist$efflen, plist$ec, plist$count, 41392, detail = TRUE)

## RNASeqQuant GD
gdest <- Adam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'Softmax'), list())
gdest <- Adam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'Softplus'), list())
gdest <- Adam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'ISRU'), list(alpha = 0.1))


gdest <- Adagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'Softmax'), list())
gdest <- Adagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'Softplus'), list())
gdest <- Adagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'ISRU'), list(alpha = 0.1))

gdest <- RMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'Softmax'), list())
gdest <- RMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'Softplus'), list())
gdest <- RMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'ISRU'), list(alpha = 0.1))

gdest <- NRMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'Softmax'), list())
gdest <- NRMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'Softplus'), list())
gdest <- NRMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1024, 0.1, list(af = 'ISRU'), list(alpha = 0.1))


## nice test
## Adam mini-batch
gdest <- AdaMax(plist$efflen, plist$ec, plist$count, length(plist$efflen), 500, 1024, 0.1, list(af = 'Softmax'), list())
gdest <- AdamW(plist$efflen, plist$ec, plist$count, 1/w, length(plist$efflen), 500, 1024, 0.01, list(af = 'Softmax'), list())
gdest <- AMSGrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 500, 1024, 0.1, list(af = 'Softmax'), list())
gdest <- Adam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 600, 1024, 0.25, list(af = 'Softmax'), list())
gdest <- NAdam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 600, 1024, 0.25, list(af = 'Softmax'), list())
gdest <- NAdagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 500, 1024, 0.1, list(af = 'Softmax'), list())

## NRMSProp full batch
emest <- EM(plist$efflen, plist$ec, plist$count, length(plist$efflen), detail = TRUE)
gdest <- NRMSProp(plist$efflen, plist$ec, plist$count, length(plist$efflen), 1067, 36580, 0.005, TRUE, list(af = 'Softmax'), list())
gdest <- NRMSPropW(plist$efflen, plist$ec, plist$count, 1/w, length(plist$efflen), 600, 36580, 0.005, list(af = 'Softmax'), list())
gdest <- NAdagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 1067, 36580, 0.7, FALSE, list(af = 'Softmax'), list())

## plot
tibble(iter = c(1 : length(emest$ll), 1 : length(gdest$ll)),
       ll = -c(emest$ll, gdest$ll),
       method = c(rep('EM', length(emest$ll)), rep('GD', length(gdest$ll)))) %>%
  slice(300 : length(emest$ll), (300 + length(emest$ll)) : (length(emest$ll) + length(emest$ll))) %>%
  ggplot(aes(x = iter, y = ll, group = method)) +
  geom_point(aes(color = method), size = 0.01)


## merge res
mergeres <- cbind(kallistoest, emest$counts, gdest$counts)
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
