\name{imp.cat}
\alias{imp.cat}
\title{
Impute missing categorical data
}
\description{
Performs single random imputation of missing values in a categorical
dataset under a user-supplied value of the underlying cell
probabilities.
}
\usage{
imp.cat(s, theta)
}
\arguments{
\item{s}{
summary list of an incomplete categorical dataset created by the
function \code{prelim.cat}.
}
\item{theta}{
parameter value under which the missing data are to be imputed.
This is an array of cell probabilities of dimension \code{s$d} whose
elements sum to one, such as produced by \code{em.cat}, \code{ecm.cat},
\code{da.cat}, \code{mda.cat} or \code{dabipf}.
}}
\value{
If the original incomplete dataset was in ungrouped format
(\code{s$grouped=F}), then a matrix like \code{s$x} except that all \code{NA}s have
been filled in.


If the original dataset was grouped, then a list with the following
components:



\item{x}{
Matrix of levels for categorical variables


}
\item{counts}{
vector of length \code{nrow(x)} containing frequencies or counts
corresponding to the levels in \code{x}.
}}
\details{
Missing data are drawn independently for each observational unit from
their conditional predictive distribution given the observed data and
\code{theta}. 
}
\note{
IMPORTANT: The random number generator seed must be set by the
function \code{rngseed} at least once in the current session before this
function can be used. 
}
\seealso{
\code{\link{prelim.cat}}, \code{\link{rngseed}}, \code{\link{em.cat}}, \code{\link{da.cat}}, \code{\link{mda.cat}}, \code{\link{ecm.cat}},
\code{\link{dabipf}} 
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
