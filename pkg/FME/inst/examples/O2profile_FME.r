
##################################################################
# A simple model of oxygen, diffusing along a spatial gradient,
# with imposed upper and lower boundary concentration
# oxygen is consumed at maximal fixed rate, monod limitation
##################################################################

par(mfrow=c(2,2))
require(FME)


##===============================================================##
##===============================================================##
##                           the Model                           ##
##===============================================================##
##===============================================================##


O2fun <- function(t,O2,pars)
{
with (as.list(pars),{

  Flux <- -diff(c(upO2,O2,lowO2))/dX
  dO2  <- -diff(Flux)/dx-cons*O2/(O2+ks)

  return(list(dO2,UpFlux = Flux[1],LowFlux = Flux[n+1]))
})
}

# The model parameters
pars <- c(upO2=360,  # concentration at upper boundary, mmolO2/m3
          lowO2=10,  # concentration at lower boundary, mmolO2/m3
          cons=80,   # consumption rate, mmolO2/m3/day
          ks=1)      # O2 half-saturation ct, mmolO2/m3

# The model grid
n  <- 100                       # nr grid points
dx <- 0.05   #cm
dX <- c(dx/2,rep(dx,n-1),dx/2)  # dispersion distances; half the grid size near boundaries
X  <- seq(dx/2,len=n,by=dx)     # distance from upper interface at middle of box

# Solve the steady-state conditions of the model
ox <- steady.band(y=runif(n),func=O2fun,parms=pars,nspec=1)

# Plot results
plot(ox$y,X,ylim=rev(range(X)),xlab="mmol/m3",
     main="Oxygen", ylab="depth, cm",type="l",lwd=2)

##===============================================================##
##===============================================================##
##       Global sensitivity analysis : Sensitivity ranges        ##
##===============================================================##
##===============================================================##

# 1. Sensitivity parameter range: consumption rate
pRange <- data.frame (min=60,max=100)
rownames(pRange) <- "cons"

# 2. Calculate sensitivity ranges for O2
#    model is solved 100 times, uniform parameter distribution (default)

Sens <- sensRange(parms=pars,y=runif(n),func=O2fun,nspec=1,num=100,
                  solver="steady.band",parRange=pRange,time=0)$summ

# same, now with normal distribution of consumption (mean = 80, variance=100)
Sens2 <- sensRange(parms=pars,y=runif(n),func=O2fun,nspec=1,dist="norm",time=0,
           num=100,solver="steady.band",parMean=c(cons=80),parCovar=100)$summ

# 3. Plot results
O2 <- Sens2[1:100,]    # first 100 columns are O2
o2range<-range(cbind(O2$Min,O2$Max))

plot(0,ylim=rev(range(X)),xlim=o2range,xlab= "O2",
     ylab="depth, cm",main="Sensitivity O2 model",type="n")

polygon(c(O2$Min,rev(O2$Max)),c(X,rev(X)),
        col=grey(0.9),border=NA)
polygon(c(O2$Mean-O2$Sd,rev(O2$Mean+O2$Sd)),
         c(X,rev(X)),col=grey(0.8),border=NA)
lines(O2$Mean,X,lwd=2)
legend("bottomright",fill=c(grey(0.9),grey(0.8)),
       legend=c("Min-Max","Mean+-sd"),bty="n")
legend("right",lty=1,lwd=2,legend="Mean",bty="n")


##===============================================================##
##===============================================================##
##       Local sensitivity analysis : sensitivity functions      ##
##===============================================================##
##===============================================================##

# Sensitivity functions
O2sens <- sensFun(func=O2fun,y=runif(n),parms=pars,
                  nspec=1,solver="steady.band")

# univariate sensitivity
format(O2sens$model,digits=2)

# bivariate sensitivity
panel.cor <- function(x, y)
             text(x=mean(range(x)),y=mean(range(y)),
             labels=format(cor(x,y),digits=2))

pairs(O2sens$fun[,-(1:2)],upper.panel=panel.cor)
mtext(outer=TRUE,side=3,line=-1.5,
      "Sensitivity functions",cex=1.5)

cor(O2sens$fun[,-(1:2)])

# multivariate sensitivity
Coll <- collin(O2sens$fun[,-(1:2)])
Coll
plot(Coll$N,Coll$collinearity,xlab="N",ylab="collinearity",
     main="identifiability",log="y")

abline(h=20,col="red")

##===============================================================##
##===============================================================##
##         Fitting the model to the data - using nlminb          ##
##===============================================================##
##===============================================================##

# The data
O2dat <- data.frame(x=seq(0.1,3.5,by=0.1),
    y = c(279,262,246,230,217,203,189,175,162,150,138,127,116,
          106,96,87,78,70,62,55,47,41,35,30,25,20,16,13,10,8,5,3,2,1,0))
O2depth <- cbind(name="O2",O2dat)        # oxygen versus depth
O2flux  <- c(UpFlux=170,LowFlux=0)       # measured fluxes

# 1. Objective function to minimise; all parameters are fitted

Objective <- function (x)
{
 pars[]<- x

 # Solve the steady-state conditions of the model
 ox <- steady.band(y=runif(n),func=O2fun,parms=pars,nspec=1)

 # Model cost: first the oxygen profile
 modO2 <- data.frame(x=X,O2=ox$y)
 Cost  <- modCost(obs=O2depth,model=modO2,x="x",y="y")

 # then the fluxes
 modFl <- c(UpFlux=ox$UpFlux,LowFlux=ox$LowFlux)
 Cost <- modCost(obs=O2flux,model=modFl,x=NULL,cost=Cost)

 return(Cost$model)
}

# 2. nlminb finds the minimum; parameters constrained to be > 0
print(system.time(Fit<-nlminb(start=c(360,10,80,1),
                  obj=Objective,lower=c(0,0,0,0))))

Fit

##===============================================================##
##===============================================================##
## Fit the model to data using the Levenberg-Marquardt algorithm ##
##===============================================================##
##===============================================================##

# 1. Define the model residuals

Residual <- function(xx)     # X have to be positive -> the log-transformed X's are fitted
{
  pars[]<-exp(xx)                    # parameter values
 # Solve the steady-state conditions of the model
 ox <- steady.band(y=runif(n),func=O2fun,parms=pars,nspec=1)
 # Model cost:
 modO2 <- data.frame(x=X,O2=ox$y)
 Cost  <- modCost(obs=O2depth,model=modO2,x="x",y="y")

 modFl <- c(UpFlux=ox$UpFlux,LowFlux=ox$LowFlux)
 Cost <- modCost(obs=O2flux,model=modFl,x=NULL,cost=Cost)

 return(Cost$residual$res)        # return residuals between model and data
}

# 2. nls.lm fits the model to the data - note the transformation to ensure positivity
# This algorithm also requires better initial conditions...
FitMrq <- nls.lm(par=log(c(360,10,80,1)),fn=Residual)

# 3. the best parameter values...
(Bestpar <- exp(FitMrq$par))

pars[]<-Bestpar                   # parameter values
# Solve the steady-state conditions of the model
ox <- steady.band(y=runif(n),func=O2fun,parms=pars,nspec=1)

plot(O2depth$y,O2depth$x,ylim=rev(range(O2depth$x)),pch=18,
     main="Oxygen-fitted", xlab="mmol/m3",ylab="depth, cm")
lines(ox$y,X)

# Model cost:
modO2 <- data.frame(x=X,O2=ox$y)
Cost  <- modCost(obs=O2depth,model=modO2,x="x",y="y",scaleVar=TRUE)

modFl <- c(UpFlux=ox$UpFlux,LowFlux=ox$LowFlux)
Cost <- modCost(obs=O2flux,model=modFl,x=NULL,cost=Cost,scaleVar=TRUE)

Cost

# Plot residuals
plot(Cost$residual$x, Cost$residual$res,xlab="depth",ylab="",main="residual")



