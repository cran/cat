\name{da.cat}
\alias{da.cat}
\title{
Data Augmentation algorithm for incomplete categorical data
}
\description{
Markov-Chain Monte Carlo method for simulating draws from the
observed-data posterior distribution of underlying cell probabilities
under a saturated multinomial model. May be used in conjunction with
\code{imp.cat} to create proper multiple imputations.
}
\usage{
da.cat(s, start, prior=0.5, steps=1, showits=FALSE)
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
\item{prior}{
optional array of hyperparameters specifying a Dirichlet prior
distribution. The default is the Jeffreys prior (all hyperparameters =
supplied with hyperparameters set to \code{NA} for those cells.
}
\item{steps}{
number of data augmentation steps to be taken. Each step consists of an
imputation or I-step followed by a posterior or P-step. 
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
of \code{theta} (I-step), and then a new value of \code{theta} is drawn from its
Dirichlet posterior distribution given the complete data (P-step).
After a suitable number of steps are taken, the resulting value of the
parameter may be regarded as a random draw from its observed-data
posterior distribution.


When the pattern of observed data is close to a monotone pattern, then
\code{mda.cat} is preferred because it will tend to converge more quickly.
}
\note{
IMPORTANT: The random number generator seed must be set at least once
by the function \code{rngseed} before this function can be used.
}
\references{
Schafer (1996) \emph{Analysis of Incomplete Multivariate Data,}
  Chapman \& Hall, Chapter 7.
}
\seealso{
\code{\link{prelim.cat}}, \code{\link{rngseed}}, \code{\link{mda.cat}}, \code{\link{imp.cat}}.
}
\examples{
data(crimes)
x      <- crimes[,-3]
counts <- crimes[,3]
s <- prelim.cat(x,counts)        # preliminary manipulations
thetahat <- em.cat(s)            # find ML estimate under saturated model
rngseed(7817)                    # set random number generator seed
theta <- da.cat(s,thetahat,50)   # take 50 steps from MLE
ximp  <- imp.cat(s,theta)        # impute once under theta
theta <- da.cat(s,theta,50)      # take another 50 steps
ximp  <- imp.cat(s,theta)        # impute again under new theta
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
%%FT 15/07/2003 : some changes to example.
