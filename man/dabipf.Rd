\name{dabipf}
\alias{dabipf}
\title{
Data augmentation-Bayesian IPF algorithm for incomplete categorical
data
}
\description{
Markov-Chain Monte Carlo method for simulating draws from the
observed-data posterior distribution of underlying cell probabilities
under hierarchical loglinear models. May be used in conjunction with
\code{imp.cat} to create proper multiple imputations.
}
\usage{
dabipf(s, margins, start, steps=1, prior=0.5, showits=F)
}
\arguments{
\item{s}{
summary list of an incomplete categorical dataset created by the
function \code{prelim.cat}.
}
\item{margins}{
vector describing the marginal totals to be fitted. A margin
is described by the factors not summed over, and margins are separated
by zeros.  Thus c(1,2,0,2,3,0,1,3) would indicate fitting the (1,2),
(2,3), and (1,3) margins in a three-way table, i.e., the model of no
three-way association.
}
\item{start}{
starting value of the parameter. The starting value should lie in the
interior of the parameter space for the given loglinear model. If
structural zeros are present, \code{start} should contain zeros in
those positions.
}
\item{steps}{
number of complete cycles of data augmentation-Bayesian IPF to be
performed.
}
\item{prior}{
optional array of hyperparameters specifying a Dirichlet
prior distribution. The default is the Jeffreys prior (all
hyperparameters = .5). If structural zeros are present, a prior
should be supplied with hyperparameters set to \code{NA} for those cells.
}
\item{showits}{
if \code{TRUE}, reports the iterations so the user can monitor the
progress of the algorithm.
}}
\value{
array of simulated cell probabilities that satisfy the loglinear
model. If the algorithm has converged, this will be a draw from the
actual posterior distribution of the parameters.
}
\note{
The random number generator seed must be set at least once by the
function \code{rngseed} before this function can be used.


The starting value must lie in the interior of the parameter space.
Hence, caution should be used when using a maximum likelihood estimate
(e.g., from \code{ecm.cat}) as a starting value. Random zeros in a table
may produce mle's with expected cell counts of zero. This difficulty
can be overcome by using as a starting value a posterior mode
calculated by \code{ecm.cat} with prior hyperparameters greater than one.
}
\seealso{
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
EB <- ifelse(B1==B2,1,0)
EI <- ifelse(I1==I2,1,0)
belt.frame <- cbind(belt.frame,EB,EI)
colnames(belt.frame)
a <- xtabs(Freq ~ D + S + B2 + I2 + EB + EI,
                data=belt.frame)
m <- list(c(1,2,3,4),c(3,4,5,6),c(1,5),
           c(1,6),c(2,6))
b <- loglin(a,margin=m)                  # fits (DSB2I2)B2I2EBEI)(DEB)(DEI)(SEI)
                                         # in Schafer's p. 304

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
#   Prior hyperparameter is set at 0.5 as in Shchafer's p. 329
#
imputations <- vector("list",10)

for (i in 1:10) {
cat("Doing imputation ",i,"\n")
  theta <- dabipf(s,m,theta,prior=0.5,   # toy chain; for comparison with
                   steps=25)             # results in Schafer's book the next
                                         # statement should be run,
                                         # rather than this one.
  \dontrun{theta <- dabipf(s,m,theta,prior=0.5,steps=2500)}			   
  imputations[[i]] <- imp.cat(s,theta)
}
  

		   
detach(belt.frame)
#
#  Example 2   (reproduces analysis performed in Schafer's p. 327.)
#
#  Caveat! I try to reproduce what has been done in that page, but although
#  the general appearance of the boxplots generated below is quite similar to
#  that of  Schafer's Fig. 8.4 (p. 327), the VALUES of the log odds do not
#  quite fall in line with those reported by said author. It doesn't look like
#  the difference can be traced to decimal vs. natural logs. On the other hand,
#  Fig. 8.4 refers to log odds, while the text near the end of page 327 gives
#  1.74 and 1.50 as the means of the *odds* (not log odds). FT, 22.7.2003.
#
#
data(older)                              # reading data
x      <- older[,1:6]                    # preliminary manipulations
counts <- older[,7]
s <- prelim.cat(x,counts)                
colnames(x)                              # names of columns
rngseed(1234)
m <- c(1,2,3,4,5,0,1,2,3,5,6,0,4,3)      # model (ASPMG)(ASPMD)(GD) in
                                         # Schafer's p. 327
                                         # do analysis with different priors
theta   <- ecm.cat(s,m,prior=1.5)        # Strong pull to uniform table
                                         # for initial estimates
prob1   <- dabipf(s,m,theta,steps=100,   # Burn-in period 
                  prior=0.1)
prob2   <- dabipf(s,m,theta,steps=100,   # Id. with second prior
                  prior=1.5)

lodds   <- matrix(0,5000,2)              # Where to store log odds ratios.

oddsr   <- function(x) {                 # Odds ratio of 2 x 2 table.
            o <- (x[1,1]*x[2,2])/
                   (x[1,2]*x[2,1])
            return(o)
            }

for(i in 1:5000) {                       # Now generate 5000 log odds
prob1  <- dabipf(s,m,prob1, prior=0.1)   
t1   <- apply(prob1,c(1,2),sum)          # Marginal GD table
                                         # Log odds ratio
lodds[i,1] <- log(oddsr(t1))
prob2  <- dabipf(s,m,prob2, prior=1.5)   # Id. with second prior
t2   <- apply(prob2,c(1,2),sum)         
lodds[i,2] <- log(oddsr(t2))
}
lodds  <- as.data.frame(lodds)
colnames(lodds) <- c("0.1","1.5")        # Similar to Schafer's Fig. 8.4.
boxplot(lodds,xlab="Prior hyperparameter")
title(main="Log odds ratio generated with DABIPF (5000 draws)")
summary(lodds)

}
\references{
Schafer (1996)  \emph{Analysis of Incomplete Multivariate Data.}
  Chapman \& Hall, Chapter 8.
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
%%FT 24/07/2003 : some changes to example.
