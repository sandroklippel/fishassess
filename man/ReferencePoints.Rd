\name{ReferencePointsBH}

\alias{Fmax}
\alias{F0.1}
\alias{F0.x}

\title{ 
Biological reference points from yield-per-recruit analysis.
}

\description{
Biological reference points for fisheries management obtained
from  yield-per-recruit analysis.
}

\usage{
Fmax(obj)
F0.1(obj)
F0.x(obj, x=1)
}

\arguments{
  \item{obj}{ object of the class \code{BH}. }
}

\details{

The yield-per-recruit analysis uses information on average
individual growth, natural mortality and vulnerability to
fishing.

\code{Fmax} is the fishing mortality for a given size at first
capture, which maximizes the average yield from each recruit
entering the fishery. It could be considered as a limit reference
point for the stock.

\code{F0.1} is the fishing mortality rate at which the slope of
the yield-per-recruit curve as a function of fishing mortality is
10\% of its value at the origin. The fishing mortality level
\code{F0.1} was proposed as a conservative target reference
point.
 
\code{F0.x} is a generalization of \code{F0.1}.

}

\value{ The value of the reference point. }
 
}

\references{

Caddy, J., 1998. A short review of precautionary reference points
an some proposals for their use in data-poor situations. Vol. 379
of FAO Fisheries Technical Paper. FAO, Rome. 

Caddy, J., Mahon, R., 1995. Reference points for fisheries
management. Vol. 374 of FAO Fisheries Technical Paper. FAO, Rome. 

Gulland, J., 1968. The concept of marginal yield from exploited
fish stocks. J. Cons. CIEM 32, 256-261. 
}

\note{
\code{Fmax} searches for the maximum between 0 and 100 of fishing
mortalities. \code{F0.1} is searched between zero and
\code{Fmax}. The slope is obtained by first derivative.
}

\author{ Sandro Klippel and Monica Brick Peres. }

\seealso{ \code{\link{Beverton.Holt}}, \code{\link{num.deriv}}, 
\code{\link{optimize}}, \code{\link{uniroot}} }

\examples{

BHxxx <- Beverton.Holt(F=0.5, M=1.1, Winf=286, K=0.37, 
                       t0=-0.2, Tr=0.4, Tc=1, approxim="g")
Fmax(BHxxx)
F0.1(BHxxx)
F0.x(BHxxx,x=2)

\dontrun{

## Monte Carlo approach
## This code does not run directly with example()

F <- rep(NA, 500)
for (i in 1:500){
  F[i] <- F0.1(Beverton.Holt(
  Wt=function(t){rnorm(1,mean=VBGFw(t, Winf=3000, K=0.25,
t0=-0.2),
  sd=300)}, Tr=1, Tc=5, Tmax=25, F=1, 
  M=runif(1, min=0.175, max=0.275)))
}
library(MASS)
truehist(F, xlab="F0.1", col="white")
lines(density(F), lwd=2)
}

}

\keyword{ methods }
\keyword{ optimize }
