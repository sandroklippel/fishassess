\name{NaturalMortality}

\alias{M.pauly}
\alias{M.pauly.w}
\alias{M.ralston}
\alias{M.gunderson}
\alias{M.jensen1}
\alias{M.jensen2}
\alias{M.jensen3}
\alias{M.rikhter}

\title{ Estimation of natural mortality. }

\description{ 
Estimation of natural mortality by several life-history relationships.
}

\usage{
M.pauly(Linf, K, Temp)
M.pauly.w(Winf, K, Temp)
M.ralston(K)
M.gunderson(GSI)
M.jensen1(Tmat)
M.jensen2(K)
M.jensen3(K)
M.rikhter(Tmat)
}

\arguments{
  \item{Linf}{ asymptotic length. }
  \item{Winf}{ asymptotic weight. }
  \item{K}{ growth rate. }
  \item{Temp}{ mean water temperature experienced by the stock. }
  \item{GSI}{ gonad-somatic index. }
  \item{Tmat}{ age at first maturity. }
}

\value{ Natural mortality estimate. }

\references{

Gunderson, D. R., Dygert, P. H., 1988. Reproductive effort as a predictor of natural mortality rate. 
J. Cons. int. Explor. Mer 44, 200-209.

Jensen, A., 1996. Beverton and Holt life history invariants result from optimal trade-off of 
reproduction and survival. Canadian Journal of Fisheries and Aquatic Sciences 53, 820-822.

Pauly, D., 1980. On the interrelationship between natural mortality, growth parameters and mean 
environmental temperature in 175 fish stocks. Journal du Conceil International pour l'Exploration 
de la Mer 39 (2), 175-192. 

Ralston, S., 1987. Tropical snappers and groupers. Biology and fisheries management. Westview Press, 
London, Chap. Mortality rates of snappers and groupers, pp. 375-404.

Rikhter, V., Efanov, V., 1976. On one of the approaches to estimation of natural mortality of fish 
populations. ICNAF Res. Doc. 76 (VI), 12. 

}

\author{ Sandro Klippel }

\keyword{ misc }
\keyword{ arith }
