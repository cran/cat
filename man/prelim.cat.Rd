\name{prelim.cat}
\alias{prelim.cat}
\title{
Preliminary manipulations on incomplete categorical data
}
\description{
This function performs grouping and sorting operations on 
categorical datasets with missing values. It creates a list that is
needed for input to em.cat, da.cat, imp.cat, etc.
}
\usage{
prelim.cat(x, counts, levs)
}
\arguments{
\item{x}{
categorical data matrix containing missing values. The data may be
provided either in ungrouped or grouped format. In ungrouped format,
the rows of x correspond to individual observational units, so that
nrow(x) is the total sample size. In grouped format, the rows of x
correspond to distinct covariate patterns; the frequencies are
provided through the \code{counts} argument.  In either format, the columns
correspond to variables.  The categories must be coded as consecutive
positive integers beginning with 1 (1,2,\dots), and missing values are
denoted by \code{NA}.
}
\item{counts}{
optional vector of length \code{nrow(x)} giving the frequencies corresponding
to the covariate patterns in x. The total sample size is
\code{sum(counts)}. If \code{counts} is missing, the data are assumed to be
ungrouped; this is equivalent to taking \code{counts} equal to
\code{rep(1,nrow(x))}.
}
\item{levs}{
optional vector of length \code{ncol(x)} indicating the number of levels
for each categorical variable. If missing, \code{levs[j]} is taken to be
\code{max(x[,j],na.rm=T)}.
}}
\value{
a list of seventeen components that summarize various features
of x after the data have been sorted by missingness patterns and
grouped according to the observed values. Components that might
be of interest to the user include:

\item{nmis}{
a vector of length \code{ncol(x)} containing the number of missing values
for each variable in x.
}
\item{r}{
matrix of response indicators showing the missing data patterns in x.
Dimension is (m,p) where m is number of distinct missingness patterns
in the rows of x, and p is the number of columns in x. Observed values
are indicated by 1 and missing values by 0. The row names give the
number of observations in each pattern, and the columns correspond to
the columns of x.
}
\item{d}{
vector of length \code{ncol(x)} indicating the number of levels for each
variable. The complete-data contingency table would be an array with
these dimensions. Identical to \code{levs} if \code{levs} was supplied.
}
\item{ncells}{
number of cells in the cross-classified contingency table, equal to
\code{prod(d)}.
}}
\references{
Chapters 7--8 of Schafer (1996)  \emph{Analysis of Incomplete Multivariate Data.}
  Chapman \& Hall.
}
\seealso{
  \code{\link{em.cat}},  \code{\link{ecm.cat}},  \code{\link{da.cat}}},
	 \code{\link{mda.cat}},  \code{\link{dabipf}},  \code{\link{imp.cat}}
}
\examples{
data(crimes)
crimes
s <- prelim.cat(crimes[,1:2],crimes[,3])   # preliminary manipulations
s$nmis                      # see number of missing observations per variable
s$r                         # look at missing data patterns
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
%%FT 15/07/2003 : some links added and changes to examples.





