\name{ipf}
\alias{ipf}
\title{
Iterative Proportional Fitting
}
\description{
ML estimation for hierarchical loglinear models via conventional
iterative proportional fitting (IPF).
}
\usage{
ipf(table, margins, start, eps=0.0001, maxits=50, showits=TRUE)
}
\arguments{
\item{table}{
contingency table (array) to be fit by a log-linear model. All
elements must be non-negative. 
}
\item{margins}{
vector describing the marginal totals to be fitted. A margin is described
by the factors not summed over, and margins are separated by zeros.
Thus c(1,2,0,2,3,0,1,3) would indicate fitting the (1,2), (2,3), and
(1,3) margins in a three-way table, i.e., the model of no three-way
association.
}
\item{start}{
starting value for IPF algorithm. The default is a uniform table.
If structural zeros appear in \code{table}, \code{start} should contain zeros
in those cells and ones elsewhere.
}
\item{eps}{
convergence criterion. This is the largest proportional change in an
expected cell count from one iteration to the next.  Any expected cell
count that drops below 1E-07 times the average cell probability
(1/number of non-structural zero cells) is set to zero during the
iterations.
}
\item{maxits}{
maximum number of iterations performed. The algorithm will stop if the
parameter still has not converged after this many iterations.
}
\item{showits}{
if \code{TRUE}, reports the iterations of IPF so the user can monitor the
progress of the algorithm.
}}
\value{
array like \code{table}, but containing fitted values (expected
frequencies) under the loglinear model. 
}
\section{DETAILS}{
This function is usually used to compute ML estimates for a loglinear
model. For ML estimates, the array \code{table} should contain the observed
frequencies from a cross-classified contingency table. Because this is
the "cell-means" version of IPF, the resulting fitted values will add
up to equals \code{sum(table)}. To obtain estimated cell probabilities,
rescale the fitted values to sum to one.


This function may also be used to compute the posterior mode of the
multinomial cell probabilities under a Dirichlet prior.  For a
posterior mode, set the elements of \code{table} to (observed frequencies +
Dirichlet hyperparameters - 1). Then, after running IPF, rescale the
fitted values to sum to one.
}
\note{
This function is essentially the same as the old S function \code{loglin}, but
results are computed to double precision.  See \code{help(loglin)} for more
details.
}
\references{
Agresti, A. (1990) Categorical Data Analysis. J. Wiley & Sons, New
York.

Schafer (1996) \emph{Analysis of Incomplete Multivariate Data.}
  Chapman & Hall, Chapter 8.
}

\seealso{
  \code{\link{ecm.cat}}, \code{\link{bipf}} 
}
\examples{
data(HairEyeColor)                     # load data
m=c(1,2,0,1,3,0,2,3)                   # no three-way interaction
fit1 <- ipf(HairEyeColor,margins=m,
            showits=TRUE)              # fit model
X2 <- sum((HairEyeColor-fit1)^2/fit1)  # Pearson chi square statistic
df <- prod(dim(HairEyeColor)-1)        # Degrees of freedom for this example
1 - pchisq(X2,df)                      # p-value for fit1
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
%%FT 15/07/2003 : some links added and changes to examples.
