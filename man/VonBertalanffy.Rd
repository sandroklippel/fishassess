\name{VonBertalanffy}
\alias{VBGFw}
\alias{VBGFl}
\alias{VBGFinv}

\title{ Von Bertalanffy growth function }

\description{
Weight and length based Von Bertalanffy growth function
and its inverse.
}

\usage{
VBGFw(t, Winf, K, t0, b=3)
VBGFl(t, Linf, K, t0)
VBGFinv(l, Linf, K, t0)
}

\arguments{
  \item{t}{ age. }
  \item{l}{ length. }  
  \item{Winf}{ asymptotic weight. }
  \item{Linf}{ asymptotic length. }
  \item{K}{ growth rate. }
  \item{t0}{ age at zero weight/length. }
  \item{b}{ alometric coefficient. }
}

\value{
Weight/length at age (t) or age at length.
}

\references{
Sparre, P., Venema, S., 1998. Introduction to tropical fish stock 
assessment - Part 1: Manual. Vol. 306/1 of FAO Fisheries Technical 
Paper. FAO, Rome.

Von Bertalanffy, L., 1938. A quantitative theory of organic growth 
(inquiries on growth laws. II). Human Biology 10 (2), 181-213. 
}

\author{ Sandro Klippel }

\seealso{ \code{\link{phi}}, \code{\link{phi.prime}} }

\keyword{ models }