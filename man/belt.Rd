\name{belt}
\alias{belt}
\alias{belt.frame}
\title{Data on driver injury and seat belt use}
\description{Data on driver injury and seat belt use.
}
\usage{data(belt)}
\format{
  The data frame \code{belt.frame} contains the following columns:
  \describe{
    \item{I}{Injury to driver (I1=Reported by police, I2=Follow up}
    \item{B}{Belt use (B1=Reported by police, B2=Follow up}
    \item{D}{Damage to vehicle (high, low)}
    \item{S}{Sex: Male or Female}
    \item{Freq}{Count}  }
}
\note{A matrix \code{belt} with similarly named columns exists that can be input
  directly to functions which do not admit data frames. Both the data
  frame and matrix include all complete and incomplete cases, from the
  police reports and follow up study.}
\source{
 Schafer (1996)  \emph{Analysis of Incomplete Multivariate Data.}
 Chapman \& Hall, Section 7.4.3, which cites

 Hochberg, Y. (1977) On the use of double sampling schemes in analyzing
 categorical data with misclassification errors, \emph{JASA}, vol. 71,
 p. 914-921.
}
\keyword{datasets}
