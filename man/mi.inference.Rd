\name{mi.inference}
\alias{mi.inference}
\title{
Multiple imputation inference
}
\description{
Combines estimates and standard errors from m complete-data analyses
performed on m imputed datasets to produce a single inference.  Uses
the technique described by Rubin (1987) for multiple
imputation inference for a scalar estimand.
}
\usage{
mi.inference(est, std.err, confidence=0.95)
}
\arguments{
\item{est}{
a list of $m$ (at least 2) vectors representing estimates (e.g.,
vectors of estimated regression coefficients) from complete-data
analyses performed on $m$ imputed datasets.
}
\item{std.err}{
a list of $m$ vectors containing standard errors from the
complete-data analyses corresponding to the estimates in \code{est}.
}
\item{confidence}{
desired coverage of interval estimates.
}}
\value{
a list with the following components, each of which is a vector of the
same length as the components of \code{est} and \code{std.err}:

\item{est}{
the average of the complete-data estimates.
}
\item{std.err}{
standard errors incorporating both the between and the
within-imputation uncertainty (the square root of the "total
variance").
}
\item{df}{
degrees of freedom associated with the t reference distribution used
for interval estimates.
}
\item{signif}{
P-values for the two-tailed hypothesis tests that the estimated
quantities are equal to zero.
}
\item{lower}{
lower limits of the (100*confidence)\% interval estimates.
}
\item{upper}{
upper limits of the (100*confidence)\% interval estimates.
}
\item{r}{
estimated relative increases in variance due to nonresponse.
}
\item{fminf}{
estimated fractions of missing information.
}}
\section{METHOD}{
Uses the method described on pp. 76-77 of Rubin (1987) for combining
the complete-data estimates from $m$ imputed datasets
for a scalar estimand. Significance levels and interval estimates are
approximately valid for each one-dimensional estimand, not for all of
them jointly.
}
\seealso{
  \code{\link{dabipf}},  \code{\link{imp.cat}}
}
\examples{
#
#  Example 1   Based on Schafer's p. 329 and ss. This is a toy version,
#              using a much shorter length of chain than required. To 
#              generate results comparable with those in the book, edit
#              the \dontrun{ } line below and comment the previous one.
#
data(belt)
attach(belt.frame)

oddsr   <- function(x) {                 # Odds ratio of 2 x 2 table.
            o <- (x[1,1]*x[2,2])/
                   (x[1,2]*x[2,1])
            o.sd <- sqrt(1/x[1,1] +      # large sample S.D. (Fienberg,
                      1/x[1,2] +         # p. 18)
                      1/x[2,1] +
                      1/x[2,2])
            return(list(o=o,sd=o.sd))
            }

colns <- colnames(belt.frame)

a <- xtabs(Freq ~ D + S + B2 + I2 + B1 + I1,
                data=belt.frame)
m <- list(c(1,2,5,6),c(1,2,3,4),c(3,4,5,6),
           c(1,3,5),c(1,4,6),c(2,4,6))
b <- loglin(a,margin=m)                  # fits (DSB1I1)(DSB2I2)(B2I2B1I1)(DB1B2)
                                         #  (DI1I2)(SI1I2) in Schafer's p. 329
s <- prelim.cat(x=belt[,-7],counts=belt[,7])
m <- c(1,2,5,6,0,1,2,3,4,0,3,4,5,6,0,1,3,5,0,1,4,6,0,2,4,6)
theta <- ecm.cat(s,margins=m,            # excruciantingly slow; needs 2558
                   maxits=5000)          # iterations.
rngseed(1234)
#
#   Now ten multiple imputations of the missing variables B2, I2 are
#   generated, by running a chain and taking every 2500th observation.
#   Prior hyperparameter is set at 0.5 as in Schafer's p. 329
#
est <- std.error <- vector("list",10)

for (i in 1:10) {
cat("Doing imputation ",i,"\n")
  theta <- dabipf(s,m,theta,prior=0.5,   # toy chain; for comparison with
                   steps=25)             # results in Schafer's book the next
                                         # statement should be run,
                                         # rather than this one.
  \dontrun{theta <- dabipf(s,m,theta,prior=0.5,steps=2500)}			   
  imp<- imp.cat(s,theta)
  imp.frame <- cbind(imp$x,imp$counts)
  colnames(imp.frame) <- colns
  a <- xtabs(Freq ~ B2 + I2,             # 2 x 2 table relating belt use
                    data=imp.frame)      # and injury
  print(a)
  odds <- oddsr(a)                       # odds ratio and std.dev.
  est[[i]] <- odds$o - 1                 # check deviations from 1 of
  std.error[[i]] <- odds$sd              # odds ratio
}
odds <- mi.inference(est,std.error)
print(odds)
detach(belt.frame)
}
\references{

  Fienberg, S.E. (1981) \emph{The Analysis of Cross-Classified Categorical
    Data}, MIT Press, Cambridge.
  
Rubin (1987) \emph{Multiple Imputation for Nonresponse in
    Surveys,} Wiley, New York, 
  
Schafer (1996)  \emph{Analysis of Incomplete Multivariate Data.}
Chapman & Hall, Chapter 8.
}
\keyword{models}

% Converted by Sd2Rd version 1.21.
%%FT 24/07/2003 : some changes to references and example added.
