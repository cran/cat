\name{em.cat}
\alias{em.cat}
\title{
EM algorithm for incomplete categorical data
}
\description{
Finds ML estimate or posterior mode of cell probabilities under
the saturated multinomial model.
}
\usage{
em.cat(s, start, prior=1, showits=TRUE, maxits=1000,
eps=0.0001)
}
\arguments{
\item{s}{
summary list of an incomplete categorical dataset produced by
the function \code{prelim.cat}.
}
\item{start}{
optional starting value of the parameter. This is an array with
dimensions \code{s$d} whose elements sum to one. The default starting value
is a uniform array (equal probabilities in all cells). If structural
zeros appear in the table, \code{start} should contain zeros in those
positions and nonzero (e.g. uniform) values elsewhere. 
}
\item{prior}{
optional vector of hyperparameters for a Dirichlet prior distribution.
The default is a uniform prior distribution (all hyperparameters = 1)
on the cell probabilities, which will result in maximum likelihood
estimation. If structural zeros appear in the table, a prior should be
supplied with \code{NA}s in those cells.
}
\item{showits}{
if \code{TRUE}, reports the iterations of EM so the user can monitor
the progress of the algorithm.
}
\item{maxits}{
maximum number of iterations performed. The algorithm will stop if the
parameter still has not converged after this many iterations.
}
\item{eps}{
convergence criterion. This is the largest proportional change in an
expected cell count from one iteration to the next.  Any expected
cell count that drops below 1E-07 times the average cell probability
(1/number of non-structural zero cells) is set to zero during the
iterations.
}}
\value{
array of dimension \code{s$d} containing the ML estimate or posterior mode,
assuming that EM  has converged by \code{maxits} iterations.
}
\note{
If zero cell counts occur in the observed-data table, the maximum
likelihood estimate may not be unique, and the algorithm may converge
to different stationary values depending on the starting value. Also,
if zero cell counts occur in the observed-data table, the ML estimate
may lie on the boundary of the parameter space. Supplying a prior with
hyperparameters greater than one will give a unique posterior mode in
the interior of the parameter space. Estimated probabilities for
structural zero cells will always be zero.
}
\references{
Schafer (1996)  \emph{Analysis of Incomplete Multivariate Data.}
  Chapman & Hall, Section 7.3.
}
\seealso{
\code{\link{prelim.cat}}, \code{\link{ecm.cat}}, \code{\link{logpost.cat}}
}
\examples{
data(crimes)
crimes
s <- prelim.cat(crimes[,1:2],crimes[,3])     # preliminary manipulations
thetahat <- em.cat(s)                        # mle under saturated model
logpost.cat(s,thetahat)                      # loglikelihood at thetahat
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
%%FT 15/07/2003 : some links added and changes to examples.
