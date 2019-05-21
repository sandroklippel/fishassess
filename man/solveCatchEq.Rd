\name{solveCatchEq}
\alias{solveCatchEq}

\title{ Solves catch equation }

\description{
Solves catch equation using Newton-Raphson iteration.
}

\usage{
solveCatchEq(catch, number2, M, accuracy=1e-04)
}

\arguments{
  \item{catch}{ catch in number in t. }
  \item{number2}{ number at sea in t+1. }
  \item{M}{ natural mortality rate. }
  \item{accuracy}{ desired accuracy. }
}

\value{
List with components
  \item{Number}{ number at sea. }
  \item{F}{ fishing mortality rate. }
  \item{accuracy}{ attained accuracy. }
}

\references{ 
Lassen, H., Medley, P., 2001. Virtual population analysis - 
A pratical manual for stock assessment. 
Vol. 400 of FAO Fisheries Technical Paper. FAO, Rome. 

Press, W., Teukolsky, S., Vetterling, W., Flannery, B., 1992. 
Numerical Recipes in FORTRAN 77: The Art of Scientific Computing, 
2nd Edition. Cambridge University Press, Cambridge.
}

\author{ Sandro Klippel }

\seealso{ \code{\link{vpa.age}} }

\examples{
solveCatchEq(25, 22.2, 0.2)
}

\keyword{ utilities }
\keyword{ iteration }
\keyword{ optimize }
