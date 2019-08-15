##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test codes~~~~~~~~~~~~~~~~~~~~~~~
library('tibble')
library('dplyr')
library('ggplot2')
library('Rcpp')
library('magrittr')
library('readr')

sourceCpp('../src/GD.cpp')
## sourceCpp('../src/WGD.cpp')
sourceCpp('../src/EM.cpp')
source('ec.R')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~small example~~~~~~~~~~~~~~~~~~~~~~~
plist <- list(ec = c('0,1,2', '1,2', '0,2', '0', '0,1'), count = rep(1, 5), efflen = rep(1, 3))

EM(plist$efflen, plist$ec, plist$count, spenum = 3) %>% .$counts

## test diff GD
GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'Momentum'), list(eta = 0.5, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'NAG'), list(eta = 0.5, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'Adagrad'), list(eta = 0.5, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'NAdagrad'), list(eta = 0.5, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'Adadelta'), list(gamma = 0.8)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'RMSProp'), list(eta = 0.1, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'Adam'), list(eta = 0.1, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'NAdam'), list(eta = 0.1, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'AdaMax'), list(eta = 0.1, decay = 0.03)) %>% .$counts

GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'AMSGrad'), list(eta = 0.1, decay = 0.03)) %>% .$counts

## test active functions
GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'Softmax', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03)) %>% .$counts
GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'SoftPlus', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03)) %>% .$counts
GD(plist$efflen, plist$ec, plist$count, spenum = 3, list(af = 'ISRU', opt = 'NRMSProp'), list(eta = 0.1, decay = 0.03)) %>% .$counts
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
gdest <- GD(plist$efflen, plist$ec, plist$count, length(plist$efflen), list(af = 'Softmax', opt = 'AdaMax'), list(eta = 0.1, decay = 0.001, assign0 = FALSE), batchsize = 1024) %>% .$counts

gdest <- AdamW(plist$efflen, plist$ec, plist$count, 1/w, length(plist$efflen), 500, 1024, 0.01, list(af = 'Softmax'), list())

gdest <- GD(plist$efflen, plist$ec, plist$count, length(plist$efflen), 600, 1024, list(af = 'Softmax', opt = 'AMSGrad'), list(eta = 0.07, decay = 0.003)) %>% .$counts

gdest <- GD(plist$efflen, plist$ec, plist$count, length(plist$efflen), 600, 1024, list(af = 'Softmax', opt = 'Adam'), list(eta = 0.1, decay = 0.03)) %>% .$counts

gdest <- GD(plist$efflen, plist$ec, plist$count, length(plist$efflen), 600, 1024, list(af = 'Softmax', opt = 'NAdam'), list(eta = 0.1, decay = 0.03)) %>% .$counts

gdest <- GD(plist$efflen, plist$ec, plist$count, length(plist$efflen), 600, 1024, list(af = 'Softmax', opt = 'NAdagrad'), list(eta = 0.1, decay = 0.00001)) %>% .$counts

## NRMSProp full batch
emest <- EM(plist$efflen, plist$ec, plist$count, length(plist$efflen), detail = TRUE)
gdest <- GD(plist$efflen, plist$ec, plist$count, length(plist$efflen), list(af = 'Softmax', opt = 'NRMSProp'), list(eta = 0.005, decay = 0.003, velocity = 0.95, gamma = 0.8), batchsize = 36580, details = TRUE)
gdest <- NRMSPropW(plist$efflen, plist$ec, plist$count, 1/w, length(plist$efflen), 600, 36580, 0.005, list(af = 'Softmax'), list())
gdest <- GD(plist$efflen, plist$ec, plist$count, length(plist$efflen), list(af = 'Softmax', opt = 'NAdagrad'), list(eta = 0.7, decay = 0.0001), batchsize = 36580, details = TRUE)

## plot
tibble(Iter = c(1 : length(emest$ll), 1 : length(gdest$ll)),
       NLL = -c(emest$ll, gdest$ll),
       method = c(rep('EM', length(emest$ll)), rep('Softmax', length(gdest$ll)))) %>%
  filter(Iter >= 300) %>%
  ggplot(aes(x = Iter, y = NLL, group = method)) +
  geom_point(aes(color = method), size = 0.01) +
  geom_hline(yintercept = -(emest$ll %>% .[length(.)]), linetype = 'dashed', color = 'gray')

## merge res
mergeres <- tibble(kest = kallistoest, emest = emest$counts, gdest = gdest$counts)
write_csv(mergeres, 'tmp1.csv')

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
