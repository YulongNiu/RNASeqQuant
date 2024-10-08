% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ec.R
\name{read_pseudo}
\alias{read_pseudo}
\title{Standard read in equivalence class (ec) files}
\usage{
read_pseudo(ecpath, countpath, abpath)
}
\arguments{
\item{ecpath}{The path of ec file.}

\item{countpath}{The path of count file.}

\item{abpath}{The path of abundance file.}
}
\value{
A \code{list}. 1st element is ec, 2nd element is count of ec, and the 3rd element is effective length of transcript.
}
\description{
Read in kallisto output.
}
\details{
Read in output files generated from \code{kallisto pseudo} and \code{kallisto quant} (version 0.44.0). The output files include equivalence class (ec), count, and effective length data.
}
\examples{
require('magrittr')

ecpath <- system.file('extdata', 'example_ec.ec', package = 'RNASeqQuant')
countpath <- system.file('extdata', 'example_count.tsv', package = 'RNASeqQuant')
abpath <- system.file('extdata', 'example_abundance.tsv', package = 'RNASeqQuant')

read_pseudo(ecpath, countpath, abpath)
}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
