\name{MethodsBH}
\alias{summary.BH}
\alias{plot.BH}
\alias{profile.BH}
\alias{plot.profileBH}

\title{ Methods for BH class. }

\description{
Methods for Beverton and Holt yield-per-recruit analysis.}

\usage{
summary.BH(object)
plot.BH(x, Flab="Fishing mortality", 
           YRlab="Yield/Recruit", 
           BRlab="Biomass/Recruit")
profile.BH(fitted, Fs=seq(0, round(Fmax(fitted) * 3, 0), 
                   (round(Fmax(fitted) * 3, 0)/100)))           
plot.profileBH(x, Flab="Fishing mortality", 
                  YRlab="Yield/Recruit", 
                  BRlab="Biomass/Recruit")                   
}

\arguments{
  \item{object}{ an object of the class BH which a summary is desired. }
  \item{x}{ BH object, result of Beverton.Holt. }
  \item{Flab}{ label for fishing mortality axis. }
  \item{YRlab}{ label for yield/recruit axis. }
  \item{BRlab}{ label for biomass/recruit axis. }  
  \item{fitted}{ the BH fitted model object. }
  \item{Fs}{ vector of fishing mortalities. }
}

\value{
Print or plot the results of a Beverton and Holt analysis.
}

\author{ Sandro Klippel and Monica Brick Peres }

\seealso{ \code{\link{Beverton.Holt}}, \code{\link{Fmax}} }

\examples{
BHx <- Beverton.Holt(F=0.6, M=0.2, Winf=3000, 
                     K=0.25, t0=-0.2, Tr=1, Tc=5, approxim="g")
summary(BHx)
plot(BHx)                     
}

\keyword{ methods }
\keyword{ hplot }
