% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{MatchEfflen}
\alias{MatchEfflen}
\title{Match transcript effect length}
\usage{
MatchEfflen(ec, efflenraw)
}
\arguments{
\item{ec}{A \code{std::vector<arma::uvec>} containing separated vectors, such as the output of \code{SplitEC()} in this package.}

\item{efflenraw}{A code{arma::vec} indicating the effect length of transcript.}
}
\value{
A \code{std::vector<arma::vec>} with the same length of \code{ecvec}.
}
\description{
Match transcript effect length with equivalence classes.
}
\details{
The length of \code{efflen} is equal to number of equivalence classes.
}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
\keyword{internal}
