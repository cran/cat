\name{ecm.cat}
\alias{ecm.cat}
\title{
ECM algorithm for incomplete categorical data
}
\description{
Finds ML estimate or posterior mode of cell probabilities under
a hierarchical loglinear model
}
\usage{
ecm.cat(s, margins, start, prior=1, showits=T, maxits=1000,
eps=0.0001)
}
\arguments{
\item{s}{
summary list of an incomplete categorical dataset produced by
the function \code{prelim.cat}.
}
\item{margins}{
vector describing the sufficient configurations or margins
in the desired loglinear model. A margin is described by the factors
not summed over, and margins are separated by zeros.  Thus
c(1,2,0,2,3,0,1,3) would indicate the (1,2), (2,3), and (1,3) margins
in a three-way table, i.e., the model of no three-way association.
The integers 1,2,\dots in the specified margins correspond to the
columns of the original data matrix \code{x} that was used to create \code{s}.
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
if \code{TRUE}, reports the iterations of ECM so the user can monitor
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
assuming that ECM  has converged by \code{maxits} iterations.
}
\details{
At each iteration, performs an E-step followed by a single cycle of
iterative proportional fitting.
}
\note{
If zero cell counts occur in the observed-data tables, the maximum
likelihood estimate may not be unique, and the algorithm may converge
to different stationary values depending on the starting value. Also,
if zero cell counts occur in the observed-data tables, the ML estimate
may lie on the boundary of the parameter space. Supplying a prior with
hyperparameters greater than one will give a unique posterior mode in
the interior of the parameter space. Estimated probabilities for
structural zero cells will always be zero.
}
\references{
Schafer (1996), \emph{Analysis of Incomplete Multivariate Data.}
  Chapman \& Hall, Chapter 8 

X. L. Meng and  D. B. Rubin (1991),
"IPF for contingency tables with missing data via the ECM algorithm,"
Proceedings of the Statistical Computing Section, Amer. Stat. Assoc.,
244--247.
}
\seealso{
\code{\link{prelim.cat}}, \code{\link{em.cat}}, \code{\link{logpost.cat}}
}
\examples{
data(older)                           # load data
#
#  Example 1
#
older[1:2,]                           # see partial content; older.frame also
                                      # available.
s <- prelim.cat(older[,-7],older[,7]) # preliminary manipulations
m <- c(1,2,5,6,0,3,4)                 # margins for restricted model
try(thetahat1 <- ecm.cat(s,margins=m))# will complain 
thetahat2 <- ecm.cat(s,margins=m,prior=1.1)
                                      # same model with prior information
logpost.cat(s,thetahat2)              # loglikelihood under thetahat2
#
#  Example 2   (reproduces analysis performed in Schafer's p. 327.)
#
m1 <- c(1,2,3,5,6,0,1,2,4,5,6,0,3,4)  # model (ASPMG)(ASPMD)(GD) in
                                      # Schafer's p. 327
theta1 <- ecm.cat(s,margins=m1,
                       prior=1.1)     # Prior to bring MLE away from boundary.
m2     <- c(1,2,3,5,6,0,1,2,4,5,6)    # model (ASPMG)(ASPMD)
theta2 <- ecm.cat(s,margins=m2,
                       prior=1.1)
lik1   <- logpost.cat(s,theta1)       # posterior log likelihood.
lik2   <- logpost.cat(s,theta2)       # id. for restricted model.
lrt    <- -2*(lik2-lik1)              # for testing significance of (GD)
p      <- 1 - pchisq(lrt,1)           # significance level
cat("LRT statistic for \n(ASMPG)(ASMPD) vs. (ASMPG)(ASMPD)(GD): ",lrt," with p-value = ",p)
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
%%FT 15/07/2003 : some links added and changes to examples.
