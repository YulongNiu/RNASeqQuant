% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Strsplit}
\alias{Strsplit}
\alias{SplitEC}
\title{Preprocess equivalence classes}
\usage{
Strsplit(s, delim)

SplitEC(ecraw)
}
\arguments{
\item{s}{A \code{string}.}

\item{delim}{The delimiter}

\item{ecraw}{A \code{character vector} and each element is a string with comma delimiter.}
}
\value{
\itemize{
  \item \code{Strsplit()}: A \code{arma::uvec} indicating the corresponding transcripts ID (starts from 0).
  \item \code{SplitEC}: A \code{std::vector<arma::uvec>} and each element indicates the transcripts IDs.
}
}
\description{
Split strings and equivalence classes.
}
\details{
\itemize{
  \item \code{Strsplit()}: Split a \code{string} with an user-defined delimiter.
  \item \code{SplitEC()}: Split batch of equivalence classes in \code{string} format.
}
}
\references{
\href{https://ysonggit.github.io/coding/2014/12/16/split-a-string-using-c.html}{split string in C++}
}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
\keyword{internal}
