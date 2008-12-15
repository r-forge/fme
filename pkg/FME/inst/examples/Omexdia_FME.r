##################################################################
##################################################################
######     Applications of the OMEXDIA model in simecol     ######
##################################################################
##################################################################

require(FME)
require(simecolModels)

# The data
O2data <- data.frame(x= seq(0.05,2.0,by=0.1) ,
         O2=c(267,226,188,155,125,98,76,56,40,26,15,
              8,3,0.6,0.1,0.01,0,0,0,0))

Ndata <- data.frame(x=seq(0.5,15,by=1),
NO3=c(20,10,3,1.5,0.8,0.5,0.4,0.3,0.25,0.25,0.25,0.25,0.25,0.25,0.25),
NH3=c(3,8,16,21,25,28, 29, 30, 31, 32, 30,31,33, 29, 31))

Flux <- data.frame(O2flux=600)

##===============================================================##
##===============================================================##
##                           the Model                           ##
##===============================================================##
##===============================================================##

# the model is a simecolobject OmexDiaDLL

# create an instance of the model
myOmexDia <- OmexDiaDLL()
parms(myOmexDia)["MeanFlux"]<-400

# run the model (estimates steady-state)
sol   <- out(sim(myOmexDia))$y

# plot results
Depth <- inputs(myOmexDia)$boxes$Depth
par(mfrow=c(2,2))
TOC  <- (sol$FDET+sol$SDET)*1200/10^9/2.5
plot(TOC,Depth,ylim=c(15,0),xlab="procent" ,main="TOC",
        type="l",lwd=2)
plot(sol$O2,Depth,ylim=c(15,0),xlab="mmol/m3" ,main="O2",
        type="l",lwd=2)
points(O2data$O2,O2data$x,cex=2,pch=18)
plot(sol$NO3,Depth,ylim=c(15,0),xlab="mmol/m3" ,main="NO3",
        type="l",lwd=2)
points(Ndata$NO3,Ndata$x,cex=2,pch=18)

plot(sol$NH3,Depth,ylim=c(15,0),xlim=range(c(Ndata$NH3,sol$NH3)),
         xlab="mmol/m3" ,main="NH3",
        type="l",lwd=2)
points(Ndata$NH3,Ndata$x,cex=2,pch=18)

mtext(side=3,outer=TRUE,"OmexDia-initial",line=-1.5,cex=1.5)

##===============================================================##
##===============================================================##
##       Global sensitivity analysis : Sensitivity ranges        ##
##===============================================================##
##===============================================================##
# The effects of bioturbation on profiles and oxygen fluxes...

# 1. Define Sensitivity parameter ranges
Parms <-parms(myOmexDia)

pselect<- c("dB0","mixdepth")

parRange <- data.frame(min=Parms[pselect]*0.8,max=Parms[pselect]*1.2)
rownames(parRange) <- names(Parms[pselect])

parRange

# 2. Define a function that takes as input the current parameter values (pselect)
# and generates as output part of the vertical profiles (isel) and the oxygen flux.
# sensRange expects a matrix with one row.

sFun <- function(P, isel)
{
  parms(myOmexDia)[pselect] <- P   # current parameter

  Out <- out(sim(myOmexDia))       # solve model

  # Sensitivity for O2, NO3, NH3, TOC  and O2 flux
  matrix(nr=1,c(Out$y$O2[isel],Out$y$NO3[isel],Out$y$NH3[isel],
       (Out$y$FDET[isel]+Out$y$SDET[isel])*1200/10^9/2.5,Out$O2flux))
}

# only upper 100 grid cells selected...
Npts <- 100
isel <- 1:Npts

# 3. Solve the model 100 times, select only summary statistics ($summ)
print(system.time(
Sens<-sensRange(solver=sFun,parms=Parms[pselect],dist="latin",
                parRange=parRange, isel=isel)$summ
))

head (Sens)

# 4. Plot the results...
# Split in separate profiles; names of output variables are not known...
O2  <- Sens[isel,]
NO3 <- Sens[Npts+isel,]
NH3 <- Sens[2*Npts+isel,]
TOC <- Sens[3*Npts+isel,]
O2fl<- Sens[nrow(Sens),]

# function to plot the ranges...
plotrange <- function(Conc,main="O2",xlab="conc",
                      ylab="depth,cm",Legend=FALSE)
{
 crange<-range(cbind(Conc$Min,Conc$Max))
 X <- Depth[isel]
 plot(0,ylim=rev(range(X)),xlim=crange,xlab= xlab,
     ylab=ylab,main=main,type="n")

 polygon(c(Conc$Min,rev(Conc$Max)),c(X,rev(X)),
        col=grey(0.9),border=NA)
 polygon(c(Conc$Mean-Conc$Sd,rev(Conc$Mean+Conc$Sd)),
         c(X,rev(X)),col=grey(0.8),border=NA)
 lines(Conc$Mean,X,lwd=2)
 if (Legend){
   legend("bottomright",fill=c(grey(0.9),grey(0.8)),
       legend=c("Min-Max","Mean+-sd"),bty="n")
   legend("right",lty=1,lwd=2,legend="Mean",bty="n")
 }
}

# plot for all
par(oma=c(0,0,2,0))
plotrange(O2,main="O2",Legend=TRUE)
plotrange(NO3,main="NO3")
plotrange(NH3,main="NH3")
plotrange(TOC,main="TOC")

mtext(side=3,outer=TRUE,"Sensitivity to bioturbation and mixing depth",cex=1.25)

##===============================================================##
##===============================================================##
##       Global sensitivity analysis : What-if scenarios         ##
##===============================================================##
##===============================================================##

# Effects of increasing the carbon flux on sediment-water fluxes and
# on integrated rates

# 1. Define Sensitivity parameter range
parRange <- data.frame(min=100,max=10000)
rownames(parRange) <- "MeanFlux"
parRange

# 2. Define a function that takes as input the current parameter value (P)
# and generates as output the fluxes and integrated rates,
# a matrix with one row
sFun2 <- function(P)
{
  parms(myOmexDia)["MeanFlux"] <- P   # current parameter

  Out <- out(sim(myOmexDia))          # solve model
  matrix(nr=1,c(Out$O2flux,Out$NH3flux,Out$NO3flux,Out$ODUflux,
                Out$TotMin,Out$OxicMin,Out$Denitri,Out$Nitri),
                dimnames=list(NULL,c("O2flux","NH3flux","NO3flux",
                "ODUflux","TotMin","OxicMin","Denitri","Nitri")))
}

# 3. Solve the model 100 times, parameter values regularly spaced
# keep full model output
print(system.time(
Sens<-sensRange(solver=sFun2,parms=Parms["MeanFlux"],dist="grid",
                parRange=parRange, Full=TRUE)
))
# first column is parameter value, next columns: variables
Response <- as.data.frame(Sens$sens)
head(Response)

# 4. Plot the results...

par (mfrow=c(2,2))

plot(Response$MeanFlux,y=Response$O2flux,xlab="C flux",ylab="nmol/cm2/d",
     main="O2 flux",type="l",lwd=2)
plot(Response$MeanFlux,y=Response$NH3flux,xlab="C flux",ylab="nmol/cm2/d",
     main="NH3 flux",type="l",lwd=2)
pDenit <- Response$Denitri/Response$TotMin   # part denitrified
plot(Response$MeanFlux,y=pDenit,xlab="C flux",ylab="-",
     main="fraction denitrified",type="l",lwd=2)
pOxic <- Response$OxicMin/Response$TotMin   # part denitrified
plot(Response$MeanFlux,y=pOxic,xlab="C flux", ylab="-",
     main="fraction oxic mineralisation",type="l",lwd=2)

mtext(side=3,outer=TRUE,"Sensitivity to C deposition",cex=1.25)

par(mfrow=c(1,1))

##===============================================================##
##===============================================================##
##       Local sensitivity analysis : sensitivity functions      ##
##===============================================================##
##===============================================================##

# 1. Sensitivity parameters
Parms <-parms(myOmexDia)
pselect<- c("dB0","mixdepth","MeanFlux","pFast","rFast","rSlow")

# 2. Define a function that takes as input the current parameter values (P)
# and generates as output profiles, and fluxes
# a matrix with one row
sFun3 <- function(P)
{
  parms(myOmexDia)[pselect] <- P   # current parameter

  Out <- out(sim(myOmexDia))       # solve model
  Out <- out(sim(myOmexDia))          # solve model
  # Sensitivity for O2, NO3, NH3, TOC  and O2 flux
  matrix(nr=1,c(Out$y$O2[isel],Out$y$NO3[isel],Out$y$NH3[isel],
       (Out$y$FDET[isel]+Out$y$SDET[isel])*1200/10^9/2.5,
        Out$O2flux,Out$NH3flux,Out$NO3flux,Out$ODUflux))
}

# 3. Sensitivity functions
isel<-seq(3,by=5,to=200)  # assume we have measured concentrations every 0.5 cm
Depth[isel]

Sens <- sensFun( solver=sFun3,parms=Parms[pselect])

# 4. univariate sensitivity
format(Sens$model,digits=2)

# 5. bivariate sensitivity
panel.cor <- function(x, y)
             text(x=mean(range(x)),y=mean(range(y)),
             labels=format(cor(x,y),digits=2))
pairs(Sens$fun[,-(1:2)],upper.panel=panel.cor)
mtext(outer=TRUE,side=3,
      "Sensitivity functions",cex=1.5)

cor(Sens$fun[,-(1:2)])

# 6. multivariate sensitivity
Coll <- collin(Sens$fun[,-(1:2)])
head(Coll)
tail(Coll)

plot(Coll$N,Coll$collinearity,xlab="N",ylab="collinearity",
     main="identifiability",log="y")
abline(h=20,col="red")

# 'identifiable parameter combinations' (>2 parameters)
Coll[which(Coll$N>2 & Coll$collinearity < 20),]

##===============================================================##
##===============================================================##
##         Fitting the model to the data - using nlminb          ##
##===============================================================##
##===============================================================##

# One parameter is fitted: MeanFlux

# 1. Define an objective function (to be minimised); the model cost

Objective <- function(X)
{
 parms(myOmexDia)["MeanFlux"]<-X               # set parameter values
 Out <- out(sim(myOmexDia))

 # select model output that corresponds to data
 O2flux <- Out$O2flux                          # modeled oxygen flux
 ModProf <-cbind(x=Depth,Out$y)                # vertical profiles

 # observed data are in 3 data sets.
 # compare to first set of data: oxygen flux
 Cost<- modCost(model=c(O2flux=O2flux),obs=Flux,x=NULL)

 # update model cost with second set of data, O2 profile
 Cost<- modCost(model=ModProf,obs=O2data,x="x",cost=Cost)

 # return SSR between model and data
 Cost<- modCost(model=ModProf,obs=Ndata,x="x",cost=Cost)$model

 return(Cost)
}

# 2. nlminb finds the minimum; parameters constrained to be > 0
print(system.time(Fit<-nlminb(start=100,
                  obj=Objective,lower=c(0))))
Fit

##===============================================================##
##===============================================================##
## Fit the model to data using the Levenberg-Marquardt algorithm ##
##===============================================================##
##===============================================================##

# 1. Define the model residuals
require(minpack.lm)
Residual <- function(X)     # X have to be positive -> the log-transformed X's are fitted
{
 parms(myOmexDia)["MeanFlux"]<-X                # set parameter values
 Out <- out(sim(myOmexDia))
 O2flux <- Out$O2flux                          # modeled oxygen flux
 ModProf <-cbind(x=Depth,Out$y)                # vertical profiles

 Cost<- modCost(model=c(O2flux=O2flux),obs=Flux,x=NULL)      # compare to first set of data
 Cost<- modCost(model=ModProf,obs=O2data,x="x",cost=Cost) # update with second
 res <- modCost(model=ModProf,obs=Ndata,x="x",cost=Cost)$residual$res
 res            # return residuals between model and data
}

# 2. nls.lm fits the model to the data
FitMrq <- nls.lm(par=100,fn=Residual)
FitMrq[]

# 3. run with best parameter value and show model cost
parms(myOmexDia)["MeanFlux"]<-FitMrq$par

Out <- out(sim(myOmexDia))

# select model output that corresponds to data
O2flux <- Out$O2flux                          # modeled oxygen flux
ModProf <-cbind(x=Depth,Out$y)                # vertical profiles

Cost<- modCost(model=c(O2flux=O2flux),obs=Flux,x=NULL)
Cost<- modCost(model=ModProf,obs=O2data,x="x",cost=Cost)
Cost<- modCost(model=ModProf,obs=Ndata,x="x",cCost=Cost)
Cost

# plot residuals...
plot(Cost$residual$res,main="residuals")

# 4. show best-fit
sol    <- Out$y
Depth  <- inputs(myOmexDia)$boxes$Depth
par(mfrow=c(2,2))
TOC  <- (sol$FDET+sol$SDET)*1200/10^9/2.5
plot(TOC,Depth,ylim=c(15,0),xlab="procent" ,main="TOC",
        type="l",lwd=2)
plot(sol$O2,Depth,ylim=c(15,0),xlab="mmol/m3" ,main="O2",
        type="l",lwd=2)
points(O2data$O2,O2data$x,cex=2,pch=18)
plot(sol$NO3,Depth,ylim=c(15,0),xlab="mmol/m3" ,main="NO3",
        type="l",lwd=2)
points(Ndata$NO3,Ndata$x,cex=2,pch=18)

plot(sol$NH3,Depth,ylim=c(15,0),xlim=range(c(Ndata$NH3,sol$NH3)),
         xlab="mmol/m3" ,main="NH3",
        type="l",lwd=2)
points(Ndata$NH3,Ndata$x,cex=2,pch=18)

mtext(side=3,outer=TRUE,"OmexDia-fitted",line=-1.5,cex=1.5)


