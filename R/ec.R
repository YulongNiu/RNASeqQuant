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
##' ecpath <- system.file('extdata', 'example_ec.ec', package = 'RNASeqQuant')
##' countpath <- system.file('extdata', 'example_count.tsv', package = 'RNASeqQuant')
##' abpath <- system.file('extdata', 'example_abundance.tsv', package = 'RNASeqQuant')
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
