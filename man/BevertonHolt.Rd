\name{BevertonHolt}
\alias{Beverton.Holt}

\title{ Beverton and Holt biomass and yield per recruit calculation. }

\description{
                         
Beverton and Holt biomass and yield per recruit calculation for a given set of
parameters. Default method is numerical integration, but also Gulland
approximation and relative yield concept are available.

}

\synopsis{

Beverton.Holt(F=NULL, M=NULL, Tr=NULL, Tc=NULL, Tmax=NULL, Wt=NULL, 
              Winf=NULL, K=NULL, t0=NULL, b=3, FZ=NULL, Lc=NULL, Linf=NULL, 
              MK=NULL, approxim=c("integration", "gulland", "relative"))

}

\usage{

Beverton.Holt(F, M, Tr, Tc, Tmax, Wt, approxim="i")
Beverton.Holt(F, M, Tr, Tc, Tmax, Winf, K, t0, b=3, approxim="i")
Beverton.Holt(F, M, Tr, Tc, Winf, K, t0, b=3, approxim="g")
Beverton.Holt(FZ, Lc, Linf, MK, approxim="r")

}

\arguments{
  \item{F}{ fishing mortality rate. }
  \item{M}{ natural mortality rate. }
  \item{Tr}{ age at recruitment. }
  \item{Tc}{ age at first capture. }
  \item{Tmax}{ maximum age attained. }
  \item{Wt}{ weight at age (function). }
  \item{Winf}{ asymptotic weight. }
  \item{K}{ growth rate. }
  \item{t0}{ age at zero length. }
  \item{b}{ exponent in length-weight relationship. }
  \item{FZ}{ exploitation rate, F/Z. }
  \item{Lc}{ length at first capture, L50. }
  \item{Linf}{asymptotic length. }
  \item{MK}{ M/K. }
  \item{approxim}{ method for estimate the integral. }
}


\value{

 An object of the class \code{BH} or \code{BH.rel} with the components
  \item{ YR }{ yield per recruit. }
  \item{ BR }{ biomass per recruit. }
  \item{ ... }{ input parameters. }

}
\references{ 

Beverton, R., Holt, S., 1957. On the dynamics of exploited fish populations. Vol. 2.
Sea Fish. Min. Agric. Fish Food G. B. 

Beverton, R., Holt, S., 1966. Manual of methods for fish stock assessment - Part 2.
Tables of yield functions. Vol. 38/1 de FAO Fisheries Technical Papers. FAO, Rome. 

Quinn, T., Deriso, R., 1999. Quantitative fish dynamics. Oxford University Press,
New York - Oxford. 

Sparre, P., Venema, S., 1998. Introduction to tropical fish stock assessment - Part
1: Manual. Vol. 306/1 of FAO Fisheries Technical Paper. FAO, Rome.

}

\author{ Sandro Klippel }

\seealso{ \code{\link{plot.BH}}, \code{\link{profile.BH}}, \code{\link{F0.1}},
\code{\link{Fmax}}, \code{\link{isopleths}} }

\examples{

BHx <- Beverton.Holt(F=0.6, M=0.2, Winf=3000, K=0.25, t0=-0.2, Tr=1, Tc=5, approxim="g")
BHy <- Beverton.Holt(F=0.6, M=0.2, Winf=3000, K=0.25, t0=-0.2, b=2.8, 
                     Tr=1, Tc=5, Tmax=25, approxim="i")

# It can be of observed values
TWobs <- VBGFw(1:40, Winf=3000, K=0.1, t0=-0.2, b=2.9) 
BHz <- Beverton.Holt(F=0.6, M=0.2, Wt=function(t){ TWobs[t] }, Tr=1, 
                     Tc=5, Tmax=40, approxim="i") 

\dontrun{

## Monte Carlo approach
## This code does not run directly with example()

F <- rep(NA, 500)
for (i in 1:500){
  F[i] <- F0.1(Beverton.Holt(
  Wt=function(t){rnorm(1,mean=VBGFw(t, Winf=3000, K=0.25, t0=-0.2),
  sd=300)}, Tr=1, Tc=5, Tmax=25, F=1, 
  M=runif(1, min=0.175, max=0.275)))
}
library(MASS)
truehist(F, xlab="F0.1", col="white")
lines(density(F), lwd=2)
}

}

\keyword{ methods }
\keyword{ list }
