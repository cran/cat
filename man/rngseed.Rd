\name{rngseed}
\alias{rngseed}
\title{Initialize random number generator seed}
\description{Seeds the random number generator}
\usage{rngseed(seed)}
\arguments{\item{seed}{
a positive number, preferably a large integer.
}}
\value{
\code{NULL}.
}
\note{
The random number generator seed must be set at least once
by this function before the simulation or imputation functions
in this package (\code{da.cat}, \code{imp.cat}, etc.) can be used.
}
\keyword{distribution}
% Converted by Sd2Rd version 1.21.

