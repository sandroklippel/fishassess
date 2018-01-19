\name{num.deriv}
\alias{num.deriv}

\title{ Numerical derivative by Ridders' method. }

\description{
Evaluates numerically the derivative of 
a function at a point by Ridders' method.
}

\usage{
num.deriv(fun, x, h=0.5, ...)
}

\arguments{
  \item{fun}{ function }
  \item{x}{ location of function }
  \item{h}{ initial stepsize }
  \item{\dots}{ further parameters to function }
}

\value{
  List with components
  \item{deriv}{ the derivative value. }
  \item{error}{ estimate of the error. }
}

\references{
Press, W., Teukolsky, S., Vetterling, W., Flannery, B., 1992. 
Numerical Recipes in FORTRAN 77: The Art of Scientific Computing, 
2nd Edition. Cambridge University Press, Cambridge.
}

\author{ Sandro Klippel }

\seealso{ \code{\link{F0.1}} }

\examples{

# F(x)  = x^3 / 3 - x
# F'(x) = x^2 - 1 
num.deriv(fun=function(x){x^3/3-x},-1)
(-1^2)-1

# F(x)  = sin(x)
# F'(x) = cos(x)
num.deriv(sin, 4)
cos(4)

# F(x)  = sqrt(x)
# F'(x) = 1 / (2*sqrt(x))
num.deriv(sqrt,0.1)
1/(2*sqrt(0.1))

}
\keyword{ math }
