\name{bipf}
\alias{bipf}
\title{
  Bayesian Iterative Proportional Fitting (BIPF)
}
\description{
Markov-Chain Monte Carlo method for simulating posterior draws of cell
probabilities under a hierarchical loglinear model
}
\usage{
bipf(table,margins, prior=0.5, start, steps=1, showits=FALSE)
}
\arguments{
\item{table }{
contingency table (array) to be fitted by a log-linear model. All
elements must be non-negative.
}
\item{margins }{
vector describing the marginal totals to be fitted. A margin
is described by the factors not summed over, and margins are separated
by zeros.  Thus c(1,2,0,2,3,0,1,3) would indicate fitting the (1,2),
(2,3), and (1,3) margins in a three-way table, i.e., the model of no
three-way association.
}
\item{prior }{
optional array of hyperparameters specifying a Dirichlet
prior distribution. The default is the Jeffreys prior (all
hyperparameters = .5). If structural zeros appear in \code{table}, a prior
should be supplied with hyperparameters set to \code{NA} for those cells.
}
\item{start }{
starting value for the algorithm. The default is a uniform table.  If
structural zeros appear in \code{table}, \code{start} should contain zeros in
those cells and ones elsewhere.
}
\item{steps }{
number of cycles of Bayesian IPF to be performed.
}
\item{showits}{
if \code{TRUE}, reports the iterations so the user can monitor the progress
of the algorithm.
}}
\value{
array like \code{table}, but containing simulated cell probabilities that
satisfy the loglinear model. If the algorithm has converged, this will
be a draw from the actual posterior distribution of the parameters.
}
\note{
The random number generator seed must be set at least once by the
function \code{rngseed} before this function can be used.


The starting value must lie in the interior of the parameter space.
Hence, caution should be used when using a maximum likelihood 
estimate (e.g., from \code{ipf}) as a starting value. Random zeros in a 
table may produce mle's with expected cell counts of zero, and any 
zero in a starting value is interpreted by \code{bipf} as a structural 
zero. This difficulty can be overcome by using as a starting value
calculated by \code{ipf} after adding a small positive constant (e.g.,
1/2) to each cell.
}
\seealso{
\code{\link{ipf}} and \code{\link{rngseed}}.
}
\examples{
data(HairEyeColor)                     # load data
m=c(1,2,0,1,3,0,2,3)                   # no three-way interaction
thetahat <- ipf(HairEyeColor,margins=m,
            showits=TRUE)              # fit model
thetahat <- ipf(HairEyeColor+.5,m)     # find an interior starting value 
rngseed(1234567)                       # set random generator seed
theta <- bipf(HairEyeColor,m,
              start=thetahat,prior=0.5,
              steps=50)                # take 50 steps
}
\references{
Schafer (1996) \emph{Analysis of Incomplete Multivariate Data.}
  Chapman \& Hall,  Chapter 8.
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
%%FT 15/07/2003 : some changes to example.
