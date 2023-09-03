\name{mda.cat}
\alias{mda.cat}
\title{
Monotone Data Augmentation algorithm for incomplete categorical data
}
\description{
Markov-Chain Monte Carlo method for simulating draws from the
observed-data posterior distribution of underlying cell probabilities
under a saturated multinomial model. May be used in conjunction with
\code{imp.cat} to create proper multiple imputations. Tends to converge
more quickly than \code{da.cat} when the pattern of observed data is nearly
monotone. 
}
\usage{
mda.cat(s, start, steps=1, prior=0.5, showits=FALSE)
}
\arguments{
\item{s}{
summary list of an incomplete categorical dataset created by the
function \code{prelim.cat}.
}
\item{start}{
starting value of the parameter. This is an array of cell
probabilities of dimension \code{s$d}, such as one created by \code{em.cat}.
If structural zeros appear in the table, starting values for those
cells should be zero.
}
\item{steps}{
number of data augmentation steps to be taken. Each step consists of
an imputation or I-step followed by a posterior or P-step.
}
\item{prior}{
optional vector of hyperparameters specifying a Dirichlet prior
distribution. The default is the Jeffreys prior (all hyperparameters =
supplied with hyperparameters set to \code{NA} for those cells.
}
\item{showits}{
if \code{TRUE}, reports the iterations so the user can monitor the
progress of the algorithm.
}}
\value{
an array like \code{start} containing simulated cell probabilities.
}
\details{
At each step, the missing data are randomly imputed under their
predictive distribution given the observed data and the current value
of \code{theta} (I-step) Unlike \code{da.cat}, however, not all of the missing
data are filled in, but only enough to complete a monotone pattern.
Then a new value of \code{theta} is drawn from its Dirichlet posterior
distribution given the monotone data (P-step).  After a suitable
number of steps are taken, the resulting value of the parameter may be
regarded as a random draw from its observed-data posterior
distribution.


For good performance, the variables in the original data matrix \code{x}
(which is used to create \code{s}) should be ordered according to their
rates of missingness from most observed (in the first columns) to
least observed (in the last columns).
}
\note{
IMPORTANT: The random number generator seed must be set at least once
by the function \code{rngseed} before this function can be used.
}
\references{
Schafer (1996)  \emph{Analysis of Incomplete Multivariate Data.}
  Chapman & Hall,  Chapter 7.
}
\seealso{
\code{\link{prelim.cat}}, \code{\link{rngseed}}, \code{\link{da.cat}}, \code{\link{imp.cat}}.
}
\examples{
data(older)
x      <- older[1:80,1:4]               # subset of the data with
counts <- older[1:80,7]                 # monotone pattern.
s <- prelim.cat(x,counts)               # preliminary manipulations
thetahat <- em.cat(s)                   # mle under saturated model
rngseed(7817)                           # set random generator seed
theta <- mda.cat(s,thetahat,50)         # take 50 steps from mle
ximp <- imp.cat(s,theta)                # impute under theta
theta <- mda.cat(s,theta,50)            # take another 50 steps
ximp <- imp.cat(s,theta)                # impute under new theta
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
%%FT 17/07/2003 : some changes to example.
