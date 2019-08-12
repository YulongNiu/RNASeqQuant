##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test codes~~~~~~~~~~~~~~~~~~~~~~~
library('tibble')
library('dplyr')
library('ggplot2')
library('Rcpp')
library('magrittr')

sourceCpp('../src/utilities.cpp')
sourceCpp('../src/likelihood.cpp')
sourceCpp('../src/GD.cpp')
sourceCpp('../src/WGD.cpp')
sourceCpp('../src/EM.cpp')
source('ec.R')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~small example~~~~~~~~~~~~~~~~~~~~~~~
plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))

EM(plist$efflen, plist$ec, plist$count, spenum = 3) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, list(af = 'Softmax', opt = 'Momentum'), list(eta = 0.5, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, list(af = 'Softmax', opt = 'NAG'), list(eta = 0.5, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, list(af = 'Softmax', opt = 'Adagrad'), list(eta = 0.5, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, list(af = 'Softmax', opt = 'Adam'), list(eta = 0.1, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, list(af = 'Softmax', opt = 'Adadelta'), list(eta = 0.1, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, list(af = 'Softmax', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, list(af = 'Softplus', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03)) %>% .$counts
GD(plist$efflen, plist$ec, plist$count, spenum = 3, 100, 1024, list(af = 'ISRU', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03)) %>% .$counts
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
## nice test
gdest <- AdaMax(plist$efflen, plist$ec, plist$count, length(plist$efflen), 500, 1024, 0.1, list(af = 'Softmax'), list())
gdest <- AdamW(plist$efflen, plist$ec, plist$count, 1/w, length(plist$efflen), 500, 1024, 0.01, list(af = 'Softmax'), list())
gdest <- AMSGrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 500, 1024, 0.1, list(af = 'Softmax'), list())
gdest <- GD(plist$efflen, plist$ec, plist$count, length(plist$efflen), 600, 1024, list(af = 'Softmax', opt = 'Adam'), list(eta = 0.1, decay = 0.03)) %>% .$counts
gdest <- NAdam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 600, 1024, 0.25, list(af = 'Softmax'), list())
gdest <- NAdagrad(plist$efflen, plist$ec, plist$count, length(plist$efflen), 500, 1024, 0.1, list(af = 'Softmax'), list())

## NRMSProp full batch
emest <- EM(plist$efflen, plist$ec, plist$count, length(plist$efflen), detail = TRUE)
gdest <- GD(plist$efflen, plist$ec, plist$count, length(plist$efflen), 600, 36580, list(af = 'Softmax', opt = 'NRMSProp'), list(eta = 0.005, decay = 0.003)) %>% .$counts
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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## RNASeqEM:::start_profiler("profile.out")
## tmp1 <- RNASeqEM:::LogSumExpRatio1(rnorm(1:10000))
## RNASeqEM:::stop_profiler()
