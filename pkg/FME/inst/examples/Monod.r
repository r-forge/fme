require(FME)

# 1. the observations
#---------------------
Obs <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour
plot(Obs,pch=16,cex=2,xlim=c(0,400),ylim=c(0,0.15),
     xlab="mg COD/l",ylab="1/hr")

# 2. the Monod model
#---------------------
Model <- function(p,x)   return(p[1]*x/(x+p[2]))

# 3. Fitting the model to the data
#---------------------
# define the residual function
Residuals  <- function(p) (Obs$y-Model(p,Obs$x))  #... model residuals

# fit the model to data
P      <- modFit(f=Residuals,p=c(0.1,1))

# plot best-fit model
x      <-0:375
lines(x,y=Model(P$par,x))

# summary of fit
sP    <- summary(P)
sP[]
print(sP)

# estimate of parameter covariances (to update parameters) and the model variance
Covar   <- sP$cov.scaled * 2.4^2/2
s2prior <- sP$modVariance

# set nprior = 0 to avoid updating model variance
MCMC <- modMCMC(f=Residuals,p=P$par,jump=Covar,niter=10000,
                var0=s2prior,n0=NULL,updatecov=100)

plot(MCMC,Full=TRUE)
pairs(MCMC)
cor(MCMC$pars)
cov(MCMC$pars)
sP$cov.scaled

