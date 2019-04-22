context('gradient')

####################gradient for single species###############
ecpath <- system.file('extdata', 'ath_ec.ec', package = 'RNASeqQuant')
countpath <- system.file('extdata', 'ath_count.tsv', package = 'RNASeqQuant')
abpath <- system.file('extdata', 'ath_abundance.tsv', package = 'RNASeqQuant')

plist <- read_pseudo(ecpath, countpath, abpath)
zeroidx <- plist$count > 0
plist$count %<>% `[`(zeroidx)
plist$ec %<>% `[`(zeroidx)

kallistoest <- read.table(abpath, stringsAsFactors = FALSE, header = TRUE)[, 4]

## EM
emest <- EM(plist$efflen, plist$ec, plist$count, length(plist$efflen))

## adam
## softmax
adamISRU <- Adam(plist$efflen, plist$ec, plist$count, length(plist$efflen), 200, 1000, 1)

## test1: count number
test_that('Total number of reads in EM results', {
  expect_equal(sum(plist$count), sum(emest$counts))
})

test_that('Total number of reads in GD results', {
  expect_equal(sum(plist$count), sum(adamISRU))
})
######################################################################
