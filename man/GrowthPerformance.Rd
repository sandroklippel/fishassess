\name{phi}
\alias{phi}
\alias{phi.prime}

\title{ Munro's phi }

\description{ 
Index of overall growth performance
}

\usage{
phi(K, Winf)
phi.prime(K, Linf)
}

\arguments{
  \item{K}{ growth rate }
  \item{Winf}{ asymptotic weight }
  \item{Linf}{ asymptotic length }
}

\details{
They were defined as
\deqn{\phi = log10(K) + \frac{2}{3} \times log10(Winf)}{%
      phi = log10(K) + (2/3) * log10(Winf)}
\deqn{\phi' = log10(K) + 2 \times log10(Linf)}{%
      phi-prime = log10(K) + 2*log10(Linf)}
}

\value{phi or phi-prime index}

\references{
Munro, J., Pauly, D., 1983. A simple method for comparing 
the growth of fishes and invertebrates. Fishbyte 1 (1), 5-6. 

Sparre, P., Venema, S., 1998. Introduction to tropical fish 
stock assessment - Part 1: Manual. Vol. 306/1 of FAO Fisheries 
Technical Paper. FAO, Rome.
}

\author{ Sandro Klippel and Monica Brick Peres }

\seealso{ \code{\link{VBGFl}}, \code{\link{VBGFw}} }

\keyword{ methods }
\keyword{ misc }
