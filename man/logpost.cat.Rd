\name{logpost.cat}
\alias{logpost.cat}
\title{
  Log-posterior density for incomplete categorical data
  }
\description{
Calculates the observed-data loglikelihood or log-posterior density
for incomplete categorical data under a specified value of the
underlying cell probabilities, e.g. as resulting from em.cat or
ecm.cat.
}
\usage{
logpost.cat(s, theta, prior)
}
\arguments{
\item{s}{
summary list of an incomplete categorical dataset created by the
function \code{prelim.cat}.
}
\item{theta}{
an array of cell probabilities of dimension \code{s$d}
}
\item{prior}{
optional vector of hyperparameters for a Dirichlet prior
distribution.  The default is a uniform prior distribution (all
hyperparameters = 1) on the cell probabilities, which will result
in evaluation of the loglikelihood. If structural zeros appear in the
table, a prior should be supplied with NAs in those cells and ones
(or other hyperparameters) elsewhere.
}}
\value{
the value of the observed-data loglikelihood or log-posterior density
function at \code{theta} 
}
\details{
This is the loglikelihood or log-posterior density that ignores the
missing-data mechanism.
}
\references{
Schafer (1996) \emph{Analysis of Incomplete Multivariate Data.}
  Chapman \& Hall. Section 7.3.
}
\seealso{
\code{\link{prelim.cat}}, \code{\link{em.cat}}, \code{\link{ecm.cat}}
}
\examples{
data(older)                            # load data
older[1:2,c(1:4,7)]                    # see partial content; older.frame also
                                       # available.
s <- prelim.cat(older[,1:4],older[,7]) # preliminary manipulations
m <- c(1,2,0,3,4)                      # margins for restricted model
thetahat1 <- ecm.cat(s,margins=m)      # mle 
logpost.cat(s,thetahat1)               # loglikelihood at thetahat1
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
%%FT 15/07/2003 : some links added and changes to examples.