\name{older}
\alias{older}
\alias{older.frame}
\title{Older people dataset}
\description{Data from the Protective Services Project for Older Persons
}
\usage{data(older)}
\format{
  The data frame \code{older.frame} contains the following columns:
  \describe{
    \item{M}{Mental status}
    \item{P}{ysical status}
    \item{D}{Survival status (deceased or not)}
    \item{G}{Group membership: E=experimental, C=control)}
    \item{A}{Age: Under75 and 75+}
    \item{S}{Sex: Male or Female}
    \item{Freq}{Count}  }
}
\note{A matrix \code{older} with similarley named columns exists that can be input
  directly to functions which do not admit data frames.}
\source{
 Schafer (1996)  \emph{Analysis of Incomplete Multivariate Data.}
  Chapman & Hall, Section 7.3.5.
}
\keyword{datasets}
