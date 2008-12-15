##################################################################
######    Applications of the BOD + O2 in a river model     ######
##################################################################
# Biochemical Oxygen Demand (BOD) and oxygen (O2) dynamics
# in a river - example from function steady.1D in rootSolve

require(FME)
par(mfrow=c(2,2))

# The observed data

Data<- matrix(nc=2,byrow=TRUE,data =c(
   50,   9.5,  1050,   3.5,  2050,   5.3,  3050,  10.5,
 4050,  35.7,  5050, 112.1,  6050, 180.6,  7050, 236.1,
 8050, 268.2,  9050, 284.8))
colnames(Data) <- c("x","O2")

##===============================================================##
##===============================================================##
##                           the Model                           ##
##===============================================================##
##===============================================================##

#==================#
# Model equations
#==================#
O2BOD <- function(t,state,pars)

{
with(as.list(pars),

{
  BOD <- state[1:N]
  O2  <- state[(N+1):(2*N)]

# BOD dynamics
  FluxBOD <-  v*c(BOD_0,BOD)  # fluxes due to water flow (advection)
  FluxO2  <-  v*c(O2_0,O2)

  BODrate <- r*BOD*O2/(O2+ks)  # 1-st order consumption, Monod in oxygen

#rate of change = flux gradient - consumption  + reaeration (O2)
  dBOD         <- -diff(FluxBOD)/dx  - BODrate
  dO2          <- -diff(FluxO2)/dx   - BODrate + p*(O2sat-O2)

  return(list(c(dBOD=dBOD,dO2=dO2),BODrate=BODrate))
 })
 }    # END O2BOD

#==================#
# Model parameters
#==================#
pars <- c(
v       = 1e2,       # velocity, m/day
r       = 1,         # /day, first-order decay of BOD
p       = 1,         # /day, air-sea exchange rate
ks      = 1,         # mmol O2/m3, half-saturation conc
O2sat   = 300,       # mmol/m3 saturated oxygen conc
O2_0    = 50 ,       # mmol/m3 riverine oxygen conc
BOD_0   = 1500)      # mmol/m3 riverine BOD concentration

#==================#
# Model grid
#==================#
dx      = 100       # grid size, meters
x       <- seq(dx/2,10000,by=dx)  # m, distance from river
N       <- length(x)

#==================#
# Model solution
#==================#
# initial guess
state <- c(rep(200,N),rep(200,N))

# steady-state solution
out   <- steady.1D (y=state,func=O2BOD,parms=pars, nspec=2,pos=TRUE)

# initial oxygen concentration
O2_in <- out$y[(N+1):(2*N)]

plot(x,O2_in,xlab= "Distance from river",
     ylab="Oxygen, mmol/m3",main="O2-BOD model",type="l")
     
##===============================================================##
##===============================================================##
##       Global sensitivity analysis : Sensitivity ranges        ##
##===============================================================##
##===============================================================##

# 1. Sensitivity parameter ranges
parRange=matrix(nc=2,data=c(500,600,0.1,0.5,0.5,2),byrow=TRUE)
rownames(parRange) <- c("v","r","p")

# 2. Calculate sensitivity ranges for O2 (variable (N+1):(2*N))
#    model is solved 100 times
print(system.time(
Sens<-sensRange(func=O2BOD,y=state,time=0,parms=pars,
                sensvar=(N+1):(2*N),solver="steady.1D",
                nspec=2,pos=TRUE,dist="unif",parRange=parRange,num=100)$summ
))

# 3. Plot the results...
yrange<-range(cbind(Sens$Min,Sens$Max))

plot(x,xlim=range(x),ylim=yrange,xlab= "Distance from river",
     ylab="Oxygen, mmol/m3",main="Sensitivity to v, r, p",type="n")

polygon(c(x,rev(x)),c(Sens$Min,rev(Sens$Max)),
        col=grey(0.9),border=NA)
polygon(c(x,rev(x)),c(Sens$Mean-Sens$Sd,
        rev(Sens$Mean+Sens$Sd)),col=grey(0.8),border=NA)
lines(x,Sens$Mean,lwd=2)
legend("bottomright",fill=c(grey(0.9),grey(0.8)),
       legend=c("Min-Max","Mean+-sd"),bty="n")
legend("right",lty=1,lwd=2,legend="Mean",bty="n")

##===============================================================##
##===============================================================##
##       Global sensitivity analysis : What-if scenarios         ##
##===============================================================##
##===============================================================##

# The effect of river flow on the oxygen concentration at the river outflow.

# 1. Define Sensitivity parameter range
parRange <- data.frame(min=10,max=5000)
rownames(parRange) <- "v"
parRange

sens <-  sensRange(func=O2BOD,y=rep(1,2*N),time=0,parms=pars,dist="grid",
                  sensvar=2*N,solver="steady.1D", Full=TRUE,
                  nspec=2,pos=TRUE,parRange=parRange,num=100)$sens

head(sens)
plot(sens,xlab="velocity,m/day",ylab="mmol/m3",main="oxygen at outflow")

##===============================================================##
##===============================================================##
##       Local sensitivity analysis : sensitivity functions      ##
##===============================================================##
##===============================================================##

# All parameters, all variables...
Sens <- sensFun(func=O2BOD,y=state,time=0,parms=pars,sensvar=NULL,
                solver="steady.1D",nspec=2,pos=TRUE)

# univariate sensitivity
format(Sens$model,digits=2)

# bivariate sensitivity
panel.cor <- function(x, y)
             text(x=mean(range(x)),y=mean(range(y)),
             labels=format(cor(x,y),digits=2))

pairs(Sens$fun[,-(1:2)],upper.panel=panel.cor)
mtext(outer=TRUE,side=3,line=-1.5,
      "Sensitivity functions",cex=1.5)

cor(Sens$fun[,-(1:2)])

# multivariate sensitivity
Coll <- collin(Sens$fun[,-(1:2)])
head(Coll)
tail(Coll)

plot(Coll$N,Coll$collinearity,xlab="N",ylab="collinearity",
     main="identifiability",log="y")

abline(h=20,col="red")

##===============================================================##
##===============================================================##
##         Fitting the model to the data - using nlminb          ##
##===============================================================##
##===============================================================##
# three parameters are fitted: r, p, ks

# 1. Define a model cost function
Objective <- function(X)
{
 pars[c("r","p","ks")]<-X                       # set parameter values

 out <- steady.1D (y=rep(200,2*N),func=O2BOD,   # estimate steady-state
                   parms=pars, nspec=2,pos=TRUE)

 Model <-cbind(x=x,O2=out$y[(N+1):(2*N)])       # select modeled oxygen conc

 modCost(model=Model,obs=Data,x="x")$model                # return SSR between model and data
}

# 2. nlminb finds the minimum; parameters constrained to be > 0
print(system.time(Fit<-nlminb(start=c(r=1,p=1,ks=1),
                  obj=Objective,lower=c(0,0,0))))

# 3. run the model with best-fit parameters
pars[c("r","p","ks")]<-Fit$par
out   <- steady.1D (y=state,func=O2BOD,parms=pars, nspec=2,pos=TRUE)
O2    <- out$y[(N+1):(2*N)]

# Plotting output  #
plot(Data,pch=18,cex=2,xlab= "Distance from river",
     ylab="mmol/m3",main="Oxygen",col="darkblue",ylim=c(0,300))

lines(x,O2_in,col="darkgrey")
lines(x,O2,lwd=2)
legend("bottomright",c("initial guess","fitted"),
       col=c("darkgrey","black"),lwd=c(1,2))

##===============================================================##
##===============================================================##
## Fit the model to data using the Levenberg-Marquardt algorithm ##
##===============================================================##
##===============================================================##

# 1. Define the model residuals
Residual <- function(X)     # X have to be positive -> the log-transformed X's are fitted
{
  pars[c("r","p","ks")]<-exp(X)                    # parameter values

 out   <- steady.1D (y=rep(200,2*N),func=O2BOD,    # estimate steady-state
                     parms=pars, nspec=2,pos=TRUE)
 Model <-cbind(x=x,O2=out$y[(N+1):(2*N)])          # select modeled oxygen conc

 modCost(model=Model,obs=Data,x="x")$residual$res  # return residuals between model and data
}

# 2. nls.lm fits the model to the data - note the log-transformation to ensure positivity
# This algorithm also requires better initial conditions...
FitMrq <- nls.lm(par=log(c(r=0.1,p=0.1,ks=0.1)),fn=Residual)

# 3. the best parameter values (backtransformed)...
(Bestpar <- exp(FitMrq$par))
FitMrq[]

