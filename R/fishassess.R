# FishAssess: add-on package for R system
# Copyright (C) 2001 Sandro Klippel and Monica Brick Peres
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

###################
#                 #
# VON BERTALANFFY #
#                 #
###################

##### EQ VB FOR LENGTHS ##### FORMULA
VBGFl <- function(t, Linf, K, t0){ Linf*(1-exp(-K*(t-t0))) }
##### END #####

##### EQ VB FOR WEIGHT ##### FORMULA
VBGFw <- function(t, Winf, K, t0, b=3){ Winf*(1-exp(-K*(t-t0)))^b }
##### END #####

##### EQ VB INVERTED ##### FORMULA
VBGFinv <- function(l, Linf, K, t0){ t0-(1/K)*log(1-(l/Linf)) }
##### END #####

##### NEW K AS FUNCTION OF NEW Linf #####
NewK <- function(K, Linf, Linf2){10^((log10(K)+2*log10(Linf))-(2*log10(Linf2)))}
##### END

##### GROWTH PERFORMANCE INDEX #####
phi.prime <- function(K, Linf){log10(K)+2*log10(Linf)}
##### END

##### GROWTH PERFORMANCE INDEX #####
phi <- function(K, Winf){log10(K)+(2/3)*log10(Winf)}
##### END

#####################
#                   #
# AGE BASED METHODS #
#                   #
#####################

##### VPA GULLAND ##### VPAAGE CLASS 

vpa.age <- function(agesdf, catches=NULL, M=NULL, species=NULL)
{

if (is.data.frame(agesdf))
     {  ages <- agesdf$ages
        catches <- agesdf$catches
	M <- attributes(agesdf)$M
	species <- attributes(agesdf)$species }
else {ages <- agesdf}

species <- ifelse(is.null(species),"unknown",species)

if (is.null(catches))
     {stop("catches cannot be NULL!\n")}
if (is.null(M))
     {stop("M cannot be NULL!\n")}     
if (length(ages)!=length(catches))
     {stop("ages and catches cannot have different dimensions!\n")}
    
    nclass <- length(ages)
    ages1 <- ages-1
    number <- rep(0, nclass)
    Fy <- c(rep(0, nclass-1),1)
    Fteps <- 3

cat("Converged F terminal:\n")    
while ( round(Fy[nclass],Fteps) != round(Fy[nclass-1],Fteps) )
{
    Fy[nclass] <- (Fy[nclass-1]+Fy[nclass])/2
    cat(Fy[nclass],"\n")
    for (i in seq(from=nclass, to=1, by=-1)){
       if (i==nclass){
           number[i] <- catches[i]/((Fy[i]/(Fy[i]+M))*(1-exp(-(Fy[i]+M))))
            }
       else {
           result <- solveCatchEq(catches[i],number[i+1],M)
           number[i] <- result$Number
	   Fy[i] <- result$F
            }
    }
}
    df <- data.frame(ageinf=ages1, agesup=ages, number, catch=catches, Fmortality=Fy)
    attributes(df)$title <- "Virtual Population Analysis (VPA)"
    attributes(df)$species <- species
    attributes(df)$date <- date()
    attributes(df)$M <- M
    attributes(df)$nrecruits <- number[1]
    attributes(df)$numberinsea <- sum(number)
    attributes(df)$meanF <- mean(Fy)
    attributes(df)$totalcatch <- sum(catches)
    class(df) <- "vpaage"
    df
}

##### END

##### SOLVE CATCH EQUATION - NEWTON-RAPSON METHOD #####

solveCatchEq <- function(catch, number2, M, accuracy=0.0001)
{

Number <- number2 * exp(M)                          # Only natural mortality
DeltaNumber <- catch * exp(M/2)                     # Popes's fishing mortality approximation

while ( abs(DeltaNumber) > accuracy )               # Do loop until desired accuracy
{
Number <- Number + DeltaNumber                      # Add correction to number
Z <- (log(Number)-log(number2))                     # Total mortality
fx <- ((1-(M/Z))*(Number-number2))-catch            # Calculate function
dfx <- 1-((Z-((Number-number2)/Number))*(M/(Z^2)))  # Calculate derivative
DeltaNumber <- -fx/dfx
}

Fy <- (log(Number)-log(number2)-M)
list(Number=Number, F=Fy, accuracy=DeltaNumber)

}

##### END

###########################
#                         #
# LENGTH BASED METHODS    #
#                         #
###########################

##### VPA JONES ##### VPALENGTH CLASS 

vpa.length <- function(lengthsdf, catches=NULL, Linf=NULL, K=NULL, M=NULL, a=NULL, b=NULL, species=NULL)
{

if (is.data.frame(lengthsdf))
     {  lengths <- lengthsdf$lengths
        catches <- lengthsdf$catches
	Linf <- attributes(lengthsdf)$VBGF.Linf
	K <- attributes(lengthsdf)$VBGF.K
	M <- attributes(lengthsdf)$M
	a <- attributes(lengthsdf)$LW.a
	b <- attributes(lengthsdf)$LW.b
	species <- attributes(lengthsdf)$species }
else {lengths <- lengthsdf}

a <- ifelse(is.null(a),0.00001, a)
b <- ifelse(is.null(b),3, b)
species <- ifelse(is.null(species),"",species)
VBGF.changed <- FALSE

if (is.null(catches))
     {stop("catches cannot be NULL!\n")}
if (is.null(Linf))
     {stop("Linf cannot be NULL!\n")}
if (is.null(K))
     {stop("K cannot be NULL!\n")}
if (is.null(M))
     {stop("M cannot be NULL!\n")}     
if (length(lengths)!=length(catches))
     {stop("lengths and catches cannot have different dimensions!\n")}
    
    nclass <- length(lengths)
    lclass <- lengths[nclass]-lengths[nclass-1]
    lengths2 <- c(lengths[2:nclass],NaN)
    number <- rep(0, nclass)
    FZdt <- c(rep(0, nclass-1),1)
    Fdt <- rep(0, nclass)
    Zdt <- rep(0, nclass)
    Nmeandt <- rep(NaN, nclass)
    if (lengths[nclass] > Linf){
      VBGF.changed <- TRUE
      K.old <- K
      Linf.old <- Linf
      cat("Linf is less than largest length class!\n")
      K <- NewK(K=K, Linf=Linf, Linf2=1.1*lengths[nclass])
      Linf <- 1.1*lengths[nclass]
      cat("Linf increased, and K was ajusted holding phi-prime\n")}
    dt <- (1/K)*log((Linf-lengths)/(Linf-lengths2))
    Xdt <- ((Linf-lengths)/(Linf-lengths2))^(M/(2*K))

#cat("Converged F/Z terminal:\n")
#while (round(FZdt[nclass],2) !=
#       round(((FZdt[nclass]*4)+(FZdt[nclass-1]*3)+(FZdt[nclass-2]*2)+FZdt[nclass-3])/10,2))
#{
    #FZdt[nclass] <- ((FZdt[nclass]*4)+(FZdt[nclass-1]*3)+(FZdt[nclass-2]*2)+FZdt[nclass-3])/10
    FZdt[nclass] <- 0.5
    cat(FZdt[nclass],"\n")
    for (i in seq(from=nclass, to=1, by=-1)){
       if (i==nclass){
           number[i] <- catches[i]/FZdt[i]
           Fdt[i] <- M*(FZdt[i]/(1-FZdt[i]))
           Zdt[i] <- Fdt[i]+M
                }
       else {
           number[i] <- (number[i+1]*Xdt[i]+catches[i])*Xdt[i]
           FZdt[i] <- catches[i]/(number[i]-number[i+1])
           Fdt[i] <- M*(FZdt[i]/(1-FZdt[i]))
           Zdt[i] <- Fdt[i]+M
           Nmeandt[i] <- (number[i]-number[i+1])/Zdt[i]}
            }
#}
    Wmeandt <- a*((lengths^b+lengths2^b)/2)
    Wmeandt[nclass] <- a*(lengths[nclass]+(lclass/2))^b
    Nmeandt[nclass] <- number[nclass]/Zdt[nclass]
    Bmeandt <- Wmeandt*Nmeandt
    Yield <- Wmeandt*catches
    df <- data.frame(lengths,lengths2,dt,catches,number,FZdt,Fdt,Zdt,Nmeandt,Bmeandt,Yield)
    attributes(df)$title <- "Jones' length cohort analysis"
    attributes(df)$species <- species
    attributes(df)$date <- date()
    attributes(df)$VBGF.K <- K
    attributes(df)$VBGF.Linf <- Linf
    attributes(df)$VBGF.changed <- VBGF.changed
    if (VBGF.changed){
          attributes(df)$VBGF.K.old <- K.old
          attributes(df)$VBGF.Linf.old <- Linf.old}
    attributes(df)$M <- M
    attributes(df)$LW.a <- a
    attributes(df)$LW.b <- b
    attributes(df)$nrecruits <- number[1]
    attributes(df)$numberinsea <- sum(Nmeandt)
    attributes(df)$integratedF <- sum(Fdt[1:nclass-1]*dt[1:nclass-1])
    attributes(df)$meanF <- sum(Fdt[1:nclass-1]*dt[1:nclass-1])/sum(dt[1:nclass-1])
    attributes(df)$totalbiomass <- sum(Bmeandt)
    attributes(df)$totalyield <- sum(Yield)
    class(df) <- "vpalength"
    df
}

##### END

##### SUMMARY VPA JONES ##### METHOD ##### VPALENGTH CLASS

summary.vpalength <- function(object)
{
cat("\n")
cat(attributes(object)$title,"\n")
cat("for species:",ifelse(attributes(object)$species == "",
    "unknown",attributes(object)$species),"\n")
cat("\n")
cat("Did at:",attributes(object)$date,"\n")
cat("\n")
cat("Input parameters\n")
cat("----------------\n")
cat("\n")
cat("Total yield (x10^3):",attributes(object)$totalyield/1000,"\n")
cat("Lengths from",object$lengths[1],"to",object$lengths[length(object$lengths)],"\n")
cat("\n")
cat("Natural mortality (M):",attributes(object)$M,"\n")
cat("\n")
if (attributes(object)$VBGF.changed){
cat("Von Bertalanffy growth function parameters\n")
cat("Linf =",attributes(object)$VBGF.Linf.old,"\n")
cat("K =",attributes(object)$VBGF.K.old,"\n")
cat("phi-prime =",phi.prime(Linf=attributes(object)$VBGF.Linf.old,
                 K=attributes(object)$VBGF.K.old),"\n")
cat("\n")
cat("Linf was less than largest length class!\n")
cat("Thus, Linf and K were ajusted holding phi-prime\n")
cat("Linf ajusted =",attributes(object)$VBGF.Linf,"\n")
cat("K ajusted=",attributes(object)$VBGF.K,"\n")
cat("phi-prime =",phi.prime(Linf=attributes(object)$VBGF.Linf,
                 K=attributes(object)$VBGF.K),"\n")}
else {		 
cat("Von Bertalanffy growth function parameters\n")
cat("Linf =",attributes(object)$VBGF.Linf,"\n")
cat("K =",attributes(object)$VBGF.K,"\n")
cat("phi-prime =",phi.prime(Linf=attributes(object)$VBGF.Linf,
                 K=attributes(object)$VBGF.K),"\n")}
cat("\n")
cat("Length-Weight relationship (W=a*L^b)\n")
cat("a =",attributes(object)$LW.a,"\n")
cat("b =",attributes(object)$LW.b,"\n")
cat("\n")
cat("Summary of results\n")
cat("------------------\n")
cat("\n")
cat("Total number in sea (x10^6):",attributes(object)$numberinsea/1000000,"\n")
cat("Total biomass in sea (x10^3):",attributes(object)$totalbiomass/1000,"\n")
cat("Ratio yield/biomass:",
     round(attributes(object)$totalyield/attributes(object)$totalbiomass,2),"\n")
cat("\n")
cat("Number of recruits at length",object$lengths[1],"(first class):",
    attributes(object)$nrecruits, "\n")
cat("Total fishing mortality (integrated):",
     attributes(object)$integratedF, "for",sum(object$dt[1:length(object$lengths)-1]),"year(s)\n")
cat("Mean instantaneous fishing mortality (F):",
     attributes(object)$meanF,"(per year)\n")
}

##### END

##### PRINT RESULTS OF VPA JONES ##### METHOD ##### VPALENGTH CLASS 

print.vpalength <- function(x)
{
  mat.obj <- cbind(x$lengths, x$lengths2, x$dt, x$catches,
                   x$number, x$FZdt, x$Fdt, x$Zdt, 
		   x$Nmeandt, x$Bmeandt, x$Yield)
  colnames(mat.obj) <- c("lengthinf","lengthsup","dt","catch","number","FZdt",
                         "Fdt","Zdt","Nmeandt","Bmeandt","Yield")
  rownames(mat.obj) <- rep("",length(x$lengths))
  summary.vpalength(x)
  cat("\n")
  cat("Results matrix\n")
  cat("--------------\n")
  cat("\n")
  print.matrix(mat.obj, quote=FALSE)
}

##### END

##### PLOT RESULTS OF VPA JONES ##### METHOD ##### VPALENGTH CLASS 

plot.vpalength <- function(x, xlab="Relative age", ylab="FZdt", toplab="Length", ...)
{
plot(x$FZdt ~ cumsum(x$dt), xlab=xlab, ylab=ylab, type="l", ...)
points(x$FZdt ~ cumsum(x$dt), pch=19)
axis(3, at=cumsum(x$dt), labels=c(x$lengths2[1:length(x$lengths2)-1],""))
box()
mtext(toplab, side=3, line=2.5, outer=FALSE)
}

##### END

##### LENGTH BASED THOMPSON AND BELL ##### METHOD ##### VPALENGTH CLASS 

TB.length <- function(lengths, F, factors, price=rep(0,length(lengths)), recruits,
Linf, K, M, a=0.00001, b=3, species="", compact=T)
{

df.lengths <- NULL
df.catches <- NULL
df.number <- NULL
df.FZdt <- NULL
df.Fdt <- NULL
df.Zdt <- NULL
df.Nmeandt <- NULL
df.Bmeandt <- NULL
df.Yield <- NULL
df.Value <- NULL
biomasses <- NULL
yields <- NULL
values <- NULL

tb <- function(lengths, F, price, factor, recruits, Linf, K, M, a, b)
{
nclass <- length(lengths)
lclass <- lengths[nclass]-lengths[nclass-1]
lengths2 <- c(lengths[2:nclass],NaN)
number <- rep(recruits, nclass)
catches <- rep(0, nclass)
Wmeandt <- a*((lengths^b+lengths2^b)/2)
Wmeandt[nclass] <- a*(lengths[nclass]+(lclass/2))^b
Nmeandt <- rep(NaN, nclass)
Xdt <- ((Linf-lengths)/(Linf-lengths2))^(M/(2*K))
Zdt <- (factor*F)+M
FZdt <- (factor*F)/Zdt

for (i in 2:nclass){
   number[i] <- number[i-1]*((1/Xdt[i-1]-FZdt[i-1])/(Xdt[i-1]-FZdt[i-1]))
   catches[i-1] <- (factor*F[i-1])*((number[i-1]-number[i])/Zdt[i-1])
   Nmeandt[i-1] <- (number[i-1]-number[i])/Zdt[i-1]
                   }
catches[nclass] <- number[nclass]*FZdt[nclass]
Nmeandt[nclass] <- number[nclass]/Zdt[nclass]
Bmeandt <- Nmeandt*Wmeandt
Yield <- catches*Wmeandt
Value <- Yield*price
df <- data.frame(lengths,catches,number,FZdt,Fdt=factor*F,Zdt,Nmeandt,Bmeandt,Yield,Value)
df
}

for (i in factors)
{
tmp <- tb(lengths=lengths, F=F, price=price, factor=i, recruits=recruits, Linf=Linf, K=K,
M=M, a=a, b=b)

if(is.null(df.lengths)){df.lengths <- tmp$lengths}
 else {df.lengths <- cbind(df.lengths, tmp$lengths)}

if(is.null(df.catches)){df.catches <- tmp$catches}
 else {df.catches <- cbind(df.catches, tmp$catches)}

if(is.null(df.number)){df.number <- tmp$number}
 else {df.number <- cbind(df.number, tmp$number)}

if(is.null(df.FZdt)){df.FZdt <- tmp$FZdt}
 else {df.FZdt <- cbind(df.FZdt, tmp$FZdt)}

if(is.null(df.Fdt)){df.Fdt <- tmp$Fdt}
 else {df.Fdt <- cbind(df.FZdt, tmp$Fdt)}

if(is.null(df.Zdt)){df.Zdt <- tmp$Zdt}
 else {df.Zdt <- cbind(df.Zdt, tmp$Zdt)}

if(is.null(df.Nmeandt)){df.Nmeandt <- tmp$Nmeandt}
 else {df.Nmeandt <- cbind(df.Nmeandt, tmp$Nmeandt)}

if(is.null(df.Bmeandt)){df.Bmeandt <- tmp$Bmeandt}
 else {df.Bmeandt <- cbind(df.Bmeandt, tmp$Bmeandt)}

if(is.null(df.Yield)){df.Yield <- tmp$Yield}
 else {df.Yield <- cbind(df.Yield, tmp$Yield)}
 
if(is.null(df.Value)){df.Value <- tmp$Value}
 else {df.Value <- cbind(df.Value, tmp$Value)}
 
if(is.null(biomasses)){biomasses <- sum(tmp$Bmeandt)}
 else {biomasses <- c(biomasses, sum(tmp$Bmeandt))}

if(is.null(yields)){yields <- sum(tmp$Yield)}
 else {yields <- c(yields, sum(tmp$Yield))}

if(is.null(values)){values <- sum(tmp$Value)}
 else {values <- c(values, sum(tmp$Value))}
}

df.catches <- as.data.frame(df.catches)
df.number <- as.data.frame(df.number)
df.FZdt <- as.data.frame(df.FZdt)
df.Fdt <- as.data.frame(df.Fdt)
df.Zdt <- as.data.frame(df.Zdt)
df.Nmeandt <- as.data.frame(df.Nmeandt)
df.Bmeandt <- as.data.frame(df.Bmeandt)
df.Yield <- as.data.frame(df.Yield)
df.Value <- as.data.frame(df.Value)

progressions <- data.frame(factors, biomasses, yields, values)
attributes(progressions)$title <- "Thompson and Bell analysis"
attributes(progressions)$species <- species
attributes(progressions)$date <- date()
attributes(progressions)$VBGF.K <- K
attributes(progressions)$VBGF.Linf <- Linf
attributes(progressions)$VBGF.M <- M
attributes(progressions)$LW.a <- a
attributes(progressions)$LW.b <- b

if(compact){return(progressions)}
 else
   {list(catches=df.catches, number=df.number, FZdt=df.FZdt, Fdt=df.Fdt, Zdt=df.Zdt,
Nmeandt=df.Nmeandt, Bmeandt=df.Bmeandt, Yield=df.Yield, Value=df.Value)}
}

##### END

##### MAXIMIZE YIELD FOR LENGTH BASED THOMPSON AND BELL ##### METHOD ##### VPALENGTH CLASS #####

TB.MSY <- function(lengths, F, recruits, Linf, K, M, interval=c(0,10), a=0.00001, b=3,
species="")
{

tb <- function(factor, lengths, F, recruits, Linf, K, M, a, b)
{
nclass <- length(lengths)
lclass <- lengths[nclass]-lengths[nclass-1]
lengths2 <- c(lengths[2:nclass],NaN)
number <- rep(recruits, nclass)
catches <- rep(0, nclass)
Wmeandt <- a*((lengths^b+lengths2^b)/2)
Wmeandt[nclass] <- a*(lengths[nclass]+(lclass/2))^b
Xdt <- ((Linf-lengths)/(Linf-lengths2))^(M/(2*K))
Zdt <- (factor*F)+M
FZdt <- (factor*F)/Zdt

for (i in 2:nclass){
   number[i] <- number[i-1]*((1/Xdt[i-1]-FZdt[i-1])/(Xdt[i-1]-FZdt[i-1]))
   catches[i-1] <- (factor*F[i-1])*((number[i-1]-number[i])/Zdt[i-1])
                   }
catches[nclass] <- number[nclass]*FZdt[nclass]
Yield <- catches*Wmeandt
sum(Yield)
}

opMSY <- optimize(tb, interval=interval, maximum=TRUE, lengths=lengths, F=F,
recruits=recruits, Linf=Linf, K=K, M=M, a=a, b=b)
list(species=species,MSY=opMSY$objective,MSYfactor=opMSY$maximum)
}
##### END

TB0.1 <- function(lengths, F, recruits, Linf, K, M, a=0.00001, b=3, species="", interval=c(0,10))
{
fx <- function(x, lengths, F, recruits, Linf, K, M, a, b)
{
nclass <- length(lengths)
lclass <- lengths[nclass]-lengths[nclass-1]
lengths2 <- c(lengths[2:nclass],NaN)
number <- rep(recruits, nclass)
catches <- rep(0, nclass)
Wmeandt <- a*((lengths^b+lengths2^b)/2)
Wmeandt[nclass] <- a*(lengths[nclass]+(lclass/2))^b
Xdt <- ((Linf-lengths)/(Linf-lengths2))^(M/(2*K))
Zdt <- (x*F)+M
FZdt <- (x*F)/Zdt

for (i in 2:nclass){
   number[i] <- number[i-1]*((1/Xdt[i-1]-FZdt[i-1])/(Xdt[i-1]-FZdt[i-1]))
   catches[i-1] <- (x*F[i-1])*((number[i-1]-number[i])/Zdt[i-1])
                   }
catches[nclass] <- number[nclass]*FZdt[nclass]
Yield <- catches*Wmeandt
sum(Yield)
}

DfxOrigin <- num.deriv(fun=fx, x=0, lengths=lengths, F=F, recruits=recruits, Linf=Linf, K=K, M=M, a=a, b=b)$deriv
DfxTarget <- (1/10)*DfxOrigin # 0.1

funMIN <- function(x, lengths, F, recruits, Linf, K, M, a, b, target){
Dfx <- num.deriv(fun=fx, x=x, lengths=lengths, F=F, recruits=recruits, Linf=Linf, K=K, M=M, a=a, b=b)$deriv
Dfx-target
}

result <- uniroot(f=funMIN, tol=0.0001, lengths=lengths, F=F, recruits=recruits, Linf=Linf, K=K, M=M, a=a, b=b, target=DfxTarget, interval=interval)

f01 <- result$root
yield <- fx(x=f01, lengths=lengths, F=F, recruits=recruits, Linf=Linf, K=K, M=M, a=a, b=b)

list(f01=f01, yield=yield)

#result <- optim(par=obj$F, funMIN, method="L-BFGS-B", lower=0, upper=Inf,
#                Tr=obj$Tr,Tc=obj$Tc,Tmax=obj$Tmax,M=obj$M, Wt=obj$Wt,Winf=obj$Winf,
#                K=obj$K,t0=obj$t0,b=obj$b,approxim=obj$approxim,target=DfxTarget)
#result$par[1]
}

##### END

##### MAXIMIZE ECONOMIC YIELD FOR LENGTH BASED THOMPSON AND BELL ##### METHOD
TB.MSE <- function(lengths, F, price, recruits, Linf, K, M, interval=c(0,10), a=0.00001,
b=3, species="")
{

tb <- function(factor, lengths, F, price, recruits, Linf, K, M, a, b)
{
nclass <- length(lengths)
lclass <- lengths[nclass]-lengths[nclass-1]
lengths2 <- c(lengths[2:nclass],NaN)
number <- rep(recruits, nclass)
catches <- rep(0, nclass)
Wmeandt <- a*((lengths^b+lengths2^b)/2)
Wmeandt[nclass] <- a*(lengths[nclass]+(lclass/2))^b
Xdt <- ((Linf-lengths)/(Linf-lengths2))^(M/(2*K))
Zdt <- (factor*F)+M
FZdt <- (factor*F)/Zdt

for (i in 2:nclass){
   number[i] <- number[i-1]*((1/Xdt[i-1]-FZdt[i-1])/(Xdt[i-1]-FZdt[i-1]))
   catches[i-1] <- (factor*F[i-1])*((number[i-1]-number[i])/Zdt[i-1])
                   }
catches[nclass] <- number[nclass]*FZdt[nclass]
Yield <- catches*Wmeandt
Value <- Yield*price
sum(Value)
}

opMSE <- optimize(tb, interval=interval, maximum=TRUE, lengths=lengths, F=F,
price=price, recruits=recruits, Linf=Linf, K=K, M=M, a=a, b=b)
list(species=species,MSE=opMSE$objective,MSEfactor=opMSE$maximum)
}

##### END

#######################
#                     #
# NATURAL MORTALITY   #
#                     #
#######################

##### NATURAL MORTALITY - PAULY 1980 #####
M.pauly <- function(Linf, K, Temp){
exp(-0.0152-0.279*log(Linf)+0.6543*log(K)+0.4634*log(Temp))}
##### END

##### NATURAL MORTALITY - PAULY 1980 (WEGHT) #####
M.pauly.w <- function(Winf, K, Temp){
exp(-0.4852-0.0824*log(Winf)+0.6757*log(K)+0.4627*log(Temp))}
##### END

##### NATURAL MORTALITY - RALSTON (1987)
M.ralston <- function(K){0.0666+(2.52*K)}
##### FIM

##### NATURAL MORTALITY - GUNERSON AND DYGERT (1988)
M.gunderson <- function(GSI){0.03+(1.68*GSI)}
##### FIM

##### NATURAL MORTALITY - RIKHTER AND EFANOV (1976)
M.rikhter <- function(Tmat){1.521/((Tmat^0.72)-0.155)}
##### FIM

##### NATURAL MORTALITY - JENSEN (1996)
M.jensen1 <- function(Tmat){1.65/Tmat}
##### FIM

##### NATURAL MORTALITY - JENSEN (1996)
M.jensen2 <- function(K){1.5*K}
##### FIM

##### NATURAL MORTALITY - JENSEN (1996)
M.jensen3 <- function(K){1.6*K}
##### FIM

#####################
#                   #
# TOTAL MORTALITY   #
#                   #
#####################

##### TOTAL MORTALITY - HOENIG 
Z.hoenig <- function(Tmax){exp(1.46-(1.01*log(Tmax)))}
##### FIM

##### TOTAL MORTALITY - JONES VAN ZALINGE ##### METHOD

Z.joneszalinge <- function(data, dlower=0, dupper=0, plot.it=T, ...)
{
range <- seq(from=1+dupper, to=length(data$lengths)-dlower)
K <- attributes(data)$VBGF.K
Linf <- attributes(data)$VBGF.Linf
CsumInv <- log(cumsum(data$catches[length(data$catches):1]))
Csum <- CsumInv[length(CsumInv):1]
LnC <- Csum[range]
LnL <- log(Linf-data$lengths[range])
regress <- lm(LnC ~ LnL)
Zk <- coef(regress)[2]
attr(Zk,"names") <- NULL
Z <- Zk*K
lst <- list(Z=Z, Linf=Linf, K=K, regress=regress, X=log(Linf-data$lengths), Y=Csum,
LnL=LnL, LnC=LnC)
if (plot.it){plot(lst$X, lst$Y, xlab="ln(Linf-L)", ylab="ln C(L, Linf)", type="p", ...);points(LnL,
LnC, pch=19);abline(coef(regress))}
lst
}

##### END

################
#              #
# SELECTIVITY  #
#              #
################

##### SELECTION FROM VPA ##### METHOD ##### VPALENGTH OR VPAAGE CLASS 

selection <- function(data){
l <- length(data$lengths)
lengths <- data$lengths[1:l-1]
S <- (data$Fdt[1:l-1]/data$dt[1:l-1])/max(data$Fdt[1:l-1]/data$dt[1:l-1])
list(lengths=lengths, S=S) }

##### END

##### PLOT 2 LOGISTIC FUNCTION #####

Plotlogistic2 <- function(lengths, S1, S2, D1, D2){
flogistic2 <- function(L, S1, S2, D1, D2){(1/(1+exp(S1-(S2*L))))*(1/(1+exp(D1-(D2*L))))}
L50 <- S1/S2
D50 <- D1/D2
S <- flogistic2(seq(from=min(lengths), to=max(lengths)), S1=S1, S2=S2, D1=D1,
D2=D2)
plot(c(min(lengths),max(lengths)),c(-0.1,1.1), type="n", xlab="L", ylab="S")
lines(c(L50,L50),c(0,1))
lines(c(D50,D50),c(0,1))
lines(seq(from=min(lengths), to=max(lengths)), S)
text(L50,-0.1,paste(round(L50,1)))
text(D50,-0.1,paste(round(D50,1)))
return(NULL)
}

##### END

##### LOGISTIC FUNCTION - DOME SHAPED ##### METHOD ##### SELECTION CLASS

S.f2logistic2 <- function(dat, plot.it=T, maxiter=50){
library(nls)
library(MASS)
nls.control(maxiter=maxiter)
D1ini <- mean(dat$lengths)*(log(3)/sd(dat$lengths))
D2ini <- D1ini/mean(dat$lengths)
S1ini <- (mean(dat$lengths)*(log(3)/sd(dat$lengths)))+D1ini
S2ini <- S1ini/mean(dat$lengths)
flogistic2 <<- function(L, S1, S2, D1,
D2){(1/(1+exp(S1-(S2*L))))*(1/(1+exp(D1-(D2*L))))}
regress <- nls(S ~ flogistic2(lengths, S1, S2, D1, D2), data=dat, start=list(S1=S1ini,
S2=S2ini,
D1=D1ini, D2=D2ini))
S1 <- coef(regress)[1]
attr(S1,"names") <- NULL
S2 <- coef(regress)[2]
attr(S2,"names") <- NULL
D1 <- coef(regress)[3]
attr(D1,"names") <- NULL
D2 <- coef(regress)[4]
attr(D2,"names") <- NULL
detach(package:nls)
detach(package:MASS)
L50 <- S1/S2
L75 <- (log(3)+S1)/S2
D50 <- D1/D2
D75 <- (log(3)+D1)/D2
S <- flogistic2(seq(from=min(dat$lengths), to=max(dat$lengths)), S1, S2, D1, D2)
if (plot.it){plot(c(min(dat$lengths),max(dat$lengths)),c(-0.1,1.1), type="n", xlab="L",
ylab="S");points(dat$lengths, dat$S);lines(seq(from=min(dat$lengths),
to=max(dat$lengths)),
S);lines(c(L50,L50),c(0,1));lines(c(D50,D50),c(0,1));text(L50,-0.1,paste(round(L50,1)));text(D50,-0.1,paste(round(D50,1)))}
list(S=S, lengths=dat$lengths, L50=L50, L75=L75, D50=D50, D75=D75, S1=S1,
S2=S2,
D1=D1, D2=D2)
}

##### END

##### LOGISTIC FUNCTION - DOME SHAPED ##### METHOD ##### SELECTION CLASS

S.f2logistic <- function(dat, plot.it=T, ...){
library(nls)
library(MASS)
D1ini <- mean(dat$lengths)*(log(3)/sd(dat$lengths))
D2ini <- D1ini/mean(dat$lengths)
S1ini <- (mean(dat$lengths)*(log(3)/sd(dat$lengths)))+D1ini
S2ini <- S1ini/mean(dat$lengths)
flogistic2 <<- function(L, S1, S2, D1,
D2){(1/(1+exp(S1-(S2*L))))*(1/(1+exp(D1-(D2*L))))}
regress <- nls(S ~ flogistic2(lengths, S1, S2, D1, D2), data=dat, start=list(S1=S1ini,
S2=S2ini, D1=D1ini, D2=D2ini),...)
S1infsup <- confint(regress, "S1")
attr(S1infsup,"names") <- NULL
S1inf <- S1infsup[1]
S1sup <- S1infsup[2]
S2infsup <- confint(regress, "S2")
attr(S2infsup,"names") <- NULL
S2inf <- S2infsup[1]
S2sup <- S2infsup[2]
S1 <- coef(regress)[1]
attr(S1,"names") <- NULL
S2 <- coef(regress)[2]
attr(S2,"names") <- NULL
D1infsup <- confint(regress, "D1")
attr(D1infsup,"names") <- NULL
D1inf <- D1infsup[1]
D1sup <- D1infsup[2]
D2infsup <- confint(regress, "D2")
attr(D2infsup,"names") <- NULL
D2inf <- D2infsup[1]
D2sup <- D2infsup[2]
D1 <- coef(regress)[3]
attr(D1,"names") <- NULL
D2 <- coef(regress)[4]
attr(D2,"names") <- NULL
detach(package:nls)
detach(package:MASS)
L50 <- S1/S2
L75 <- (log(3)+S1)/S2
D50 <- D1/D2
D75 <- (log(3)+D1)/D2
S <- flogistic2(seq(from=min(dat$lengths), to=max(dat$lengths)), S1, S2, D1, D2)
if (plot.it){plot(c(min(dat$lengths),max(dat$lengths)),c(-0.1,1.1), type="n", xlab="L",
ylab="S");points(dat$lengths, dat$S);lines(seq(from=min(dat$lengths),
to=max(dat$lengths)),
S);lines(c(L50,L50),c(0,1));lines(c(D50,D50),c(0,1));text(L50,-0.1,paste(round(L50,1)));t
ext(D50,-0.1,paste(round(D50,1)))}
list(S=S, lengths=dat$lengths, L50=L50, L75=L75, D50=D50, D75=D75, S1=S1,
S1inf=S1inf, S1sup=S1sup, S2=S2, S2inf=S2inf, S2sup=S2sup, D1=D1, D1inf=D1inf,
D1sup=D1sup, D2=D2, D2inf=D2inf, D2sup=D2sup)
}

##### END

##### LOGISTIC FUNCTION - `S' SELECTION ##### METHOD ##### SELECTION CLASS

S.flogistic <- function(dat, plot.it=T){
library(nls)
library(MASS)
S1ini <- 2*(mean(dat$lengths)*(log(3)/sd(dat$lengths)))
S2ini <- S1ini/mean(dat$lengths)
flogistic <<- function(L, S1, S2){1/(1+exp(S1-(S2*L)))}
regress <- nls(S ~ flogistic(lengths, S1, S2), data=dat, start=list(S1=S1ini, S2=S2ini))
S1infsup <- confint(regress, "S1")
attr(S1infsup,"names") <- NULL
S1inf <- S1infsup[1]
S1sup <- S1infsup[2]
S2infsup <- confint(regress, "S2")
attr(S2infsup,"names") <- NULL
S2inf <- S2infsup[1]
S2sup <- S2infsup[2]
S1 <- coef(regress)[1]
attr(S1,"names") <- NULL
S2 <- coef(regress)[2]
attr(S2,"names") <- NULL
detach(package:nls)
detach(package:MASS)
L50 <- S1/S2
L75 <- (log(3)+S1)/S2
S <- flogistic(seq(from=min(dat$lengths), to=max(dat$lengths)), S1, S2)
if(plot.it){plot(c(min(dat$lengths),max(dat$lengths)),c(-0.1,1.1), type="n", xlab="L",
ylab="S");points(dat$lengths, dat$S);lines(seq(from=min(dat$lengths),
to=max(dat$lengths)), S);lines(c(L50,L50),c(0,1));text(L50,-0.1,paste(round(L50,1)))}
list(S=S, lengths=dat$lengths, L50=L50, L75=L75, S1=S1, S1inf=S1inf, S1sup=S1sup,
S2=S2, S2inf=S2inf, S2sup=S2sup)
}

##### END

###################
#                 #
# BEVERTON & HOLT #
#                 #
###################

##### YIELD AND BIOMASS PER RECRUIT ##### BH CLASS #####

Beverton.Holt <- function(F=NULL,M=NULL,Tr=NULL,Tc=NULL,Tmax=NULL,Wt=NULL,Winf=NULL,K=NULL,t0=NULL,
                          b=3,FZ=NULL,Lc=NULL,Linf=NULL,MK=NULL,
                          approxim=c("integration","gulland","relative"))
{
##### RELATIVE BEVERTON & HOLT
rel.Beverton.Holt <- function(FZ, Lc, Linf, MK)
{
if (Lc > Linf | FZ > 1){
lst <- list(YR=NA, FZ=FZ, Lc=Lc, Linf=Linf, MK=MK, approxim="relative")
class(lst) <- "rel.BH"
return(lst)
}
U <- 1-(Lc/Linf)
m <- (1-FZ)/MK
YR <- FZ*(U^m)*(1-((3*U)/(1+m))+((3*U^2)/(1+2*m))-((U^3)/(1+3*m)))
lst <- list(YR=YR, FZ=FZ, Lc=Lc, Linf=Linf, MK=MK, approxim="relative")
class(lst) <- "rel.BH"
lst
}

# test inputs
if (substring(approxim,1,1) == "r"){
if (any(c(is.null(FZ),is.null(Lc),is.null(Linf),is.null(MK))))
   { stop("FZ, Lc, Linf, MK parameters are required for relative yield approximation") }
return(rel.Beverton.Holt(FZ=FZ, Lc=Lc, Linf=Linf, MK=MK))
}

if ((substring(approxim,1,1) == "g" | substring(approxim,1,1) == "i") 
     & any(c(is.null(F),is.null(M),is.null(Tr),is.null(Tc))))
   { stop("F,M,Tr,Tc parameters are required for integration or Gulland approximation") }
if (substring(approxim,1,1) == "g" & any(c(is.null(t0),is.null(Winf),is.null(K))))
   { stop("VBGF parameters (Winf, K, t0) are required for Gulland approximation") }
if (substring(approxim,1,1) == "i" & is.null(Tmax)){ stop("Tmax is required for integration") }
if (substring(approxim,1,1) == "i" & !is.function(Wt) & any(c(is.null(t0),is.null(Winf),is.null(K))))
   { stop("Wt function or VBGF parameters (Winf, K, t0) are required for integration") }
if (all(c("i","g","r") != substring(approxim,1,1))){
    stop("approxim must be integration, gulland, relative")}

# gulland approximation
if (substring(approxim,1,1) == "g")
{
# cat("using Gulland approximation\n")
S <- exp(-K*(Tc-t0))
Z <- F+M
BR <- Winf*exp(-M*(Tc-Tr))*((1/Z)-((3*S)/(Z+K))+((3*(S^2))/(Z+(2*K)))-(S^3/(Z+(3*K))))
YR <- BR*F
} else 
{
if (substring(approxim,1,1) == "i" & is.function(Wt)){ 
# cat("integrate with function Wt\n")
integrand <- function(t,Tr,Tc,M,Wt){ exp((-M*(Tc-Tr))-((M+F)*(t-Tc)))*Wt(t) } 
result <- integrate(integrand, lower=Tc, upper=Tmax, stop.on.error=FALSE, Tr=Tr, Tc=Tc, M=M, Wt=Wt) 
} else {if (substring(approxim,1,1) == "i" & !is.function(Wt)){ 
       #cat("integrate with VBGF\n")
       integrand <- function(t,Tr,Tc,M,Winf,K,t0,b)
                            { exp((-M*(Tc-Tr))-((M+F)*(t-Tc)))*VBGFw(t,Winf,K,t0,b) } 
       result <- integrate(integrand, lower=Tc, subdivisions=1000, upper=Tmax, stop.on.error=FALSE, Tr=Tr, 
                           Tc=Tc, M=M, Winf=Winf, K=K, t0=t0, b=b)
       }}
if (result$message != "OK"){cat(result$message,"\n"); stop("Abort !")}
BR <- result$value
YR <- BR*F 
}

lst <- list(YR=YR, BR=BR, Tr=Tr, Tc=Tc, Tmax=Tmax, M=M, F=F, Wt=Wt,
            Winf=Winf, K=K, t0=t0, b=b, approxim=approxim[1])
class(lst) <- "BH"
lst
}

##### END

##### YIELD/RECRUIT OVER FISHING MORTALITY (PROFILE) ##### profileBH CLASS #####

profile.BH <- function(fitted, Fs=seq(0,round(Fmax(fitted)*3,0),(round(Fmax(fitted)*3,0)/100))){
if (class(fitted) != "BH"){stop("Argument must be BH class")}
F <- rep(0,length(Fs))
YR <- rep(0,length(Fs))
BR <- rep(0,length(Fs))
for (i in seq(1,length(Fs)))
{ 
result <- Beverton.Holt(F=Fs[i],Tr=fitted$Tr,Tc=fitted$Tc,Tmax=fitted$Tmax,M=fitted$M,Wt=fitted$Wt,
                        Winf=fitted$Winf,K=fitted$K,t0=fitted$t0,b=fitted$b,approxim=fitted$approxim) 
F[i] <- Fs[i]
YR[i] <- result$YR
BR[i] <- result$BR
}
lst <- list(F=F, YR=YR, BR=BR, Tr=fitted$Tr,Tc=fitted$Tc,Tmax=fitted$Tmax,M=fitted$M,Wt=fitted$Wt,
            Winf=fitted$Winf,K=fitted$K,t0=fitted$t0,b=fitted$b,approxim=fitted$approxim)
class(lst) <- "profileBH"
lst
}

##### END

##### PLOT PROFILE FOR BEVERTON & HOLT WITH FMAX E F0.1 ##### METHOD ##### BH CLASS #####

plot.BH <- function(x, Flab="Fishing mortality", YRlab="Yield/Recruit", BRlab="Biomass/Recruit"){ 
if (class(x) != "BH"){stop("Argument must be BH class")}
profobj <- profile.BH(x)
plot.profileBH(profobj,Flab=Flab,YRlab=YRlab,BRlab=BRlab)
P.Fmax <- Fmax(x) # FMAX
P.F01 <- F0.1(x)  # F0.1
YRmax <- Beverton.Holt(F=P.Fmax,Tr=x$Tr,Tc=x$Tc,Tmax=x$Tmax,M=x$M,Wt=x$Wt,Winf=x$Winf,
                        K=x$K,t0=x$t0,b=x$b,approxim=x$approxim)$YR/max(profobj$YR)
YR01 <- Beverton.Holt(F=P.F01,Tr=x$Tr,Tc=x$Tc,Tmax=x$Tmax,M=x$M,Wt=x$Wt,Winf=x$Winf,
                        K=x$K,t0=x$t0,b=x$b,approxim=x$approxim)$YR/max(profobj$YR)
points(x=P.Fmax, y=YRmax, pch="+", cex=1.5)
points(x=P.F01, y=YR01, pch="+", cex=1.5)
points(x=x$F, y=x$YR/max(profobj$YR), pch=19)
lines(x=c(-1,P.Fmax,P.Fmax), y=c(YRmax,YRmax,-1), lty="dashed")
lines(x=c(-1,P.F01,P.F01), y=c(YR01,YR01,-1), lty="dashed")
}

##### END

##### SUMMARY FOR BEVERTON & HOLT WITH FMAX E F0.1 ##### METHOD ##### BH CLASS #####

summary.BH <- function(object){ 
if (class(object) != "BH"){stop("Argument must be BH class")}
S.Fmax <- Fmax(object) # FMAX
S.F01 <- F0.1(object)  # F0.1
cat("\n")
cat("Beverton and Holt biomass and yield-per-recruit\n")
cat("-----------------------------------------------\n")
cat("\n")
cat("Input parameters\n")
cat("----------------\n")
cat("Method",switch(substring(object$approxim,1,1),
    i="..............: Integration", 
    g="..............: Gulland",
    "..............:Unknown"),
    "\n",sep="")
cat("Fishing mortality...:",object$F,"\n")
cat("Natural mortality...:",object$M,"\n")
cat("Age at recruitment..:",object$Tr,"\n") 
cat("Age at 1st capture..:",object$Tc,"\n") 
if (substring(object$approxim,1,1) == "i"){
cat("Maximum age attained:",object$Tmax,"\n")
  if (is.function(object$Wt)){
    cat("\n")
    cat("Function:\n"); print(object$Wt)}    
  else {
    cat("\n")
    cat("VBGF parameters (Winf, K, t0, b):", object$Winf, object$K, object$t0, object$b,"\n")}
} else {
cat("\n")
cat("VBGF parameters (Winf, K, t0):", object$Winf, object$K, object$t0,"\n")
}
cat("\n")
cat("Results\n")
cat("-------\n")
cat("Yield/Recruit:",round(object$YR,3),"\n")
cat("Biomass/Recruit:",round(object$BR,3),"\n")
cat("\n")
cat("Reference points\n")
cat("----------------\n")
cat("Fmax:",round(S.Fmax,2),"\n")
cat("F0.1:",round(S.F01,2),"\n")
}

##### END

##### PLOT PROFILE BEVERTON & HOLT ##### METHOD #####

plot.profileBH <- function(x, Flab="Fishing mortality", YRlab="Yield/Recruit", BRlab="Biomass/Recruit"){ 
if (class(x) != "profileBH"){stop("Argument must be profileBH class")}
decimal <- function(x){
d <- -9
while(10^d < x){d <- d + 1}
10^d
} 
byYR <- decimal(max(x$YR)-min(x$YR))/10
byBR <- decimal(max(x$BR)-min(x$BR))/10
scaleYR <- seq(from=floor(min(x$YR)), to=ceiling(max(x$YR)), by=byYR)
scaleBR <- seq(from=decimal(min(x$BR)), to=ceiling(max(x$BR)), by=byBR)
plot.new()
par(mar=c(4, 4, 1, 4)+0.25)
plot.window(xlim=c(min(x$F),max(x$F)), ylim=c(0,1))
lines(x=x$F, y=x$YR/max(x$YR))
lines(x=x$F, y=x$BR/max(x$BR))
axis(1)
axis(2, at=scaleYR/max(x$YR), labels=scaleYR)
axis(4, at=scaleBR/max(x$BR), labels=scaleBR)
box()
title(xlab=Flab, ylab=YRlab)
mtext(BRlab, side=4, line=3)
}


##### END

##### MAXIMIZE YPR AS FUNCTION OF F ##### METHOD ##### BH CLASS #####

Fmax <- function(obj){
if (class(obj) != "BH"){stop("Argument must be BH class")}

B.H <- function(F,Tr,Tc,Tmax,M,Wt,Winf,K,t0,b,approxim){
       Beverton.Holt(F=F,Tr=Tr,Tc=Tc,Tmax=Tmax,M=M,Wt=Wt,
                     Winf=Winf,K=K,t0=t0,b=b,approxim=approxim)$YR}

result <- optimize(f=B.H, interval=c(0,100), maximum=TRUE,
          tol=0.0001, Tr=obj$Tr, Tc=obj$Tc, Tmax=obj$Tmax, M=obj$M, Wt=obj$Wt, Winf=obj$Winf,
          K=obj$K, t0=obj$t0, b=obj$b, approxim=obj$approxim)

#result <- optim(par=obj$F, B.H, method="L-BFGS-B", lower=0, upper=Inf,
#                control = list(fnscale=-1),Tr=obj$Tr,Tc=obj$Tc,Tmax=obj$Tmax,M=obj$M,
#                Wt=obj$Wt,Winf=obj$Winf, K=obj$K,t0=obj$t0,b=obj$b,approxim=obj$approxim)
#result$par[1]
result$maximum
}

##### END

##### MAXIMIZE RELATIVE YPR AS FUNCTION OF F ##### METHOD ##### rel.BH CLASS #####

Emax <- function(obj){
if (class(obj) != "rel.BH"){stop("Argument must be rel.BH class")}
if (obj$Lc > obj$Linf | obj$FZ > 1){
return(NA)
}

B.H <- function(FZ, Lc, Linf, MK){
       Beverton.Holt(FZ=FZ, Lc=Lc, Linf=Linf, MK=MK, approxim="r")$YR}

result <- optimize(f=B.H, interval=c(0,1), maximum=TRUE,
          tol=0.0001, Lc=obj$Lc, Linf=obj$Linf, MK=obj$MK)

result$maximum
}

##### END

##### NUMERICAL DERIVATIVE ##### METHOD #####

# numerical derivative by Ridders' method 
# as described in Numerical Recipes, Press et al.

num.deriv <- function(fun, x, h=0.5, ...){

# tests

if (!is.function(fun)){stop("fun should be a function")}
if (h==0){stop("h must be nonzero")}

# initialize

TAB <- 10              # Size of Neville table
dStep <- 1.4           # decrement step
dStep2 <- dStep*dStep

a <- matrix(rep(NA,TAB*TAB),ncol=TAB,nrow=TAB)

a[1,1] <- (fun(x+h, ...)-fun(x-h, ...))/(2*h)
err <- 1.0e30

# first loop

for (i in seq(2,TAB)){
    h <- h/dStep

    # 1st derivative: f'(x) = [f(x+dx)-f(x-dx)]/2dx

    a[1,i] <- (fun(x+h, ...)-fun(x-h, ...))/(2*h) 
    fac <- dStep2
    
    # second loop
    
    for (j in seq(2,i)){
        a[j,i] <- (a[j-1,i]*fac-a[j-1,i-1])/(fac-1)
        fac <- fac*dStep2
        errt <- max(c(abs(a[j,i]-a[j-1,i]),abs(a[j,i]-a[j-1,i-1])))
        if (errt <= err){
           err <- errt
           ans <- a[j,i]
        }
    }
    if (abs(a[i,i]-a[i-1,i-1]) >= 2*err) break
}
lst <- list(deriv=ans, error=err)
return(lst)
}

##### END

##### F FOR x% MARGINAL GAIN OF YPR AT ORIGIN ##### METHOD ##### BH CLASS #####

F0.x <- function(obj, x=1){
if (class(obj) != "BH"){stop("Argument must be BH class")}

fx <- function(x,Tr,Tc,Tmax,M,Wt,Winf,K,t0,b,approxim)
{
Beverton.Holt(F=x,Tr=Tr,Tc=Tc,Tmax=Tmax,M=M,Wt=Wt,
              Winf=Winf,K=K,t0=t0,b=b,approxim=approxim)$YR
}

DfxOrigin <- num.deriv(fun=fx,x=0,Tr=obj$Tr,Tc=obj$Tc,Tmax=obj$Tmax,M=obj$M,Wt=obj$Wt,
                       Winf=obj$Winf,K=obj$K,t0=obj$t0,b=obj$b,approxim=obj$approxim)$deriv
DfxTarget <- (x/10)*DfxOrigin

funMIN <- function(x,Tr,Tc,Tmax,M,Wt,Winf,K,t0,b,approxim,target){
Dfx <- num.deriv(fun=fx,x=x,Tr=Tr,Tc=Tc,Tmax=Tmax,M=M,Wt=Wt,
                            Winf=Winf,K=K,t0=t0,b=b,approxim=approxim)$deriv
Dfx-target
}

result <- uniroot(f=funMIN,interval=c(0,Fmax(obj)),tol=0.0001,Tr=obj$Tr,Tc=obj$Tc,
                   Tmax=obj$Tmax,M=obj$M,Wt=obj$Wt,Winf=obj$Winf,K=obj$K,t0=obj$t0,
                   b=obj$b,approxim=obj$approxim,target=DfxTarget)
result$root

#result <- optim(par=obj$F, funMIN, method="L-BFGS-B", lower=0, upper=Inf,
#                Tr=obj$Tr,Tc=obj$Tc,Tmax=obj$Tmax,M=obj$M, Wt=obj$Wt,Winf=obj$Winf,
#                K=obj$K,t0=obj$t0,b=obj$b,approxim=obj$approxim,target=DfxTarget)
#result$par[1]
}

##### END

##### F FOR 10% MARGINAL GAIN OF YPR AT ORIGIN ##### METHOD ##### BH CLASS #####

F0.1 <- function(obj){
if (class(obj) != "BH"){stop("Argument must be BH class")}

dX <- .Machine$double.eps^0.25
fdX <- Beverton.Holt(F=dX,Tr=obj$Tr,Tc=obj$Tc,Tmax=obj$Tmax,M=obj$M,Wt=obj$Wt,Winf=obj$Winf,
                        K=obj$K,t0=obj$t0,b=obj$b,approxim=obj$approxim)$YR
DfX0 <- fdX/dX
DfX0.1 <- 0.1*DfX0

funBH <- function(F,Tr,Tc,Tmax,M,Wt,Winf,K,t0,b,approxim,Df){
dX <- .Machine$double.eps^0.25
FdX <- F + dX 
fX <- Beverton.Holt(F=F,Tr=Tr,Tc=Tc,Tmax=Tmax,M=M,Wt=Wt,Winf=Winf,K=K,t0=t0,b=b,approxim=approxim)$YR 
fXdX <- Beverton.Holt(F=FdX,Tr=Tr,Tc=Tc,Tmax=Tmax,M=M,Wt=Wt,Winf=Winf,K=K,t0=t0,b=b,approxim=approxim)$YR 
DfX <- (fXdX-fX)/dX
abs(Df-DfX)
}

result <- optimize(f=funBH,interval=c(0,Fmax(obj)),Tr=obj$Tr,Tc=obj$Tc,Tmax=obj$Tmax,M=obj$M,Wt=obj$Wt,Winf=obj$Winf,
                   K=obj$K,t0=obj$t0,b=obj$b,approxim=obj$approxim,Df=DfX0.1)
result$minimum
}

##### END

##### E FOR x% MARGINAL GAIN OF YPR AT ORIGIN ##### METHOD ##### rel.BH CLASS #####

E0.x <- function(obj, x=1){
if (class(obj) != "rel.BH"){stop("Argument must be rel.BH class")}
if (obj$Lc > obj$Linf | obj$FZ > 1){
return(NA)
}

D.BHrel <- function(E, Lc, Linf, MK){
(((1 - (Lc/Linf))^((1 - E)/MK)) - E * ((1 - (Lc/Linf))^((1 - 
    E)/MK) * (log((1 - (Lc/Linf))) * (1/MK)))) * (1 - ((3 * (1 - 
    (Lc/Linf)))/(1 + ((1 - E)/MK))) + ((3 * (1 - (Lc/Linf))^2)/(1 + 
    2 * ((1 - E)/MK))) - (((1 - (Lc/Linf))^3)/(1 + 3 * ((1 - 
    E)/MK)))) + E * ((1 - (Lc/Linf))^((1 - E)/MK)) * ((3 * (1 - 
    (Lc/Linf))^2) * (2 * (1/MK))/(1 + 2 * ((1 - E)/MK))^2 - (3 * 
    (1 - (Lc/Linf))) * (1/MK)/(1 + ((1 - E)/MK))^2 - ((1 - (Lc/Linf))^3) * 
    (3 * (1/MK))/(1 + 3 * ((1 - E)/MK))^2)
}

fx <- function(x, Lc, Linf, MK){
       Beverton.Holt(FZ=x, Lc=Lc, Linf=Linf, MK=MK, approxim="r")$YR}

DfxOrigin <- D.BHrel(E=0, Lc=obj$Lc, Linf=obj$Linf, MK=obj$MK)
DfxTarget <- (x/10)*DfxOrigin

funMIN <- function(x,Lc,Linf,MK,target){
Dfx <- D.BHrel(E=x,Lc=Lc,Linf=Linf,MK=MK)
#cat("Dfx",Dfx,"\n")
#cat("target",target,"\n")
abs(Dfx-target)
}

#result <- uniroot(f=funMIN,interval=c(0,1),tol=0.00001,
#                  Lc=obj$Lc, Linf=obj$Linf, MK=obj$MK, target=DfxTarget)
#result$root
result <- optimize(f=funMIN,interval=c(0,1),Lc=obj$Lc, Linf=obj$Linf, MK=obj$MK, target=DfxTarget)
result$minimum
}

##### END

##### E FOR 10% MARGINAL GAIN OF YPR AT ORIGIN ##### METHOD ##### rel.BH CLASS #####

E0.1 <- function(obj){E0.x(obj, x=1)}

##### END

##### ISOPLETH DIAGRAM - BEVERTON AND HOLT #####
isopleths <- function(x){
if (class(x) != "BH.rel"){stop("Argument must be BH.rel class")}

contour(x=as.numeric(rownames(dat$YR)), y=as.numeric(colnames(dat$YR)), 
        z=dat$YR, labcex=1, xlab="F/Z", ylab="Lc", main=paste("M/K =",x$MK))
}

##### END

