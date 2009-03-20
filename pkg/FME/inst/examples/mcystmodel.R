library("simecol")
require(FME)

## unit conversion
fem2mu <- 1e+6 * 1e-15
mu2fem <- 1e+9

## function that plots simulated and observed data
plot.mcyst <- function(ex, sc) {
  o <- out(sc)
  nf <- layout(matrix(c(1, 2), 2, 2, byrow=FALSE), c(3, 1), 1)
  opar <- par(mar=c(0,6.1,2,1), pch=16, cex=1, las=1)
  ylim_cells <- c(0, max(o$cells, na.omit(ex$cells)))
  ylim_mcyst  <- c(0, max(o$mcyst*fem2mu, na.omit(ex$mcyst)))
  plot(ex$time, ex$cells, axes=FALSE, col="forestgreen", ylab="",
    ylim=ylim_cells)
  mtext(expression(paste("Cells (",mL^{-1}, ")")),
    side=2, line=4, cex=1.2, las=0)
  axis(2)
  box()
  lines(o$time, o$cells)
  par(mar=c(4.1,6.1,0,1))
  plot(ex$time, ex$mcyst, col="navy", xlab = "Time (d)", ylab="",
    ylim=ylim_mcyst)
  mtext(expression(paste("MCYST (",mu~g~mL^{-1}, ")")),
    side=2, line=4, cex=1.2, las=0)
  lines(o$time, o$mcyst)# * fem2mu)
  par(opar)
}


mcyst2 <- new("odeModel",
  main = function (time, init, parms, ...) {
    x <- init
    with(as.list(parms),{
      mu.inst <- mu * x[1] * (1-x[1]/K)
      dx1 <-  mu.inst                    # microcystis
      dx2 <-  mu.inst * p  - dM * x[2]   # mcyst
      list(c(dx1, dx2))
    })
  },
  #init = c(
  #  cells = 0.8e6 ,     # Inoculum    (cells/L)
  #  mcyst = 2000 * 1e5  # Start-MCYST (fg/L)
  #),
  parms = c(
    mu = 0.2,         # 1/d (Lab)
    K  = 30e6,        # Cells per ml
    p  = 50,          # fg/cell
    dM = 0.04,        # 1/d (Orr & Jones)
    cells = 1.4e6,
    mcyst = 8.6 * 1e7
  ),
  times = c(from=0, to=100, by=1),
  solver = "lsoda",
  ## copy two "parameters" to init; this allows fitting initial values
  initfunc = function(obj) {
    init(obj)[c("cells", "mcyst")] <- parms(obj)[c("cells", "mcyst")]
    obj
  }
)

## Read the data...
dat <- read.table("expS108.txt", header=TRUE, sep="\t")
dat$mcyst[1] <- 0 ## replace NA

## Make the data compatible with model...
dat$mcyst <- dat$mcyst * mu2fem


## Simulate
mc2 <- mcyst2
sc  <- sim(mc2)
plot.mcyst(dat, sc)

## Parameters to be fitted...
whichpar <- c("p", "mu", "K", "dM")
#whichpar <- c("cells", "mcyst", "p", "mu", "K", "dM")
lower    <- c(cells = 1e5, mcyst = 0,   p = 0,   mu = 0.1, K = 10e6, dM = 0)
upper    <- c(cells = 5e6, mcyst = 5e8, p = 200, mu = 1,   K = 40e6, dM = 1)

#===========================
# 1. Fit, using fitOdeModel
#===========================

obstime   <- dat$time
times(sc) <- obstime

optdat <- data.frame(cells=dat$cells, mcyst=dat$mcyst )

## this works...
Oderes <- fitOdeModel(sc, whichpar=whichpar, obstime, optdat,  method="BFGS",
  debuglevel=0, lower = lower, upper = upper, control=list(trace=TRUE))

Oderes$par

# this works too...
Oderes2 <- fitOdeModel(sc, whichpar=whichpar, obstime, optdat,  method="PORT",
  debuglevel=0, lower = lower, upper = upper, control=list(trace=TRUE))

Oderes$par

## asign parms to the model object and simulate again
parms(sc)[whichpar] <- Oderes$par

times(sc) <- c(from=0, to=100, by=1)
sc        <- sim(sc)
plot.mcyst(dat, sc)

#> parms(sc)
#          mu            K            p           dM        cells        mcyst
#2.324523e-01 3.180356e+07 1.255524e+02 4.413899e-02 1.498322e+06 7.332897e+07



#===========================
# 2. Fit using other algorithms ...
#===========================
# Solver function
Solver <- function (p) {
 parms(mc2)[whichpar]<- p
 out(sim(mc2))
}

# Cost of model versus data
Cost <- function (p) {
 Run <- Solver(p)
 return ( modCost(model=Run,obs=dat,weight="std"))
}

print(system.time(Fit<-modFit(p=parms(mc2)[whichpar],
                  f=Cost,lower=lower[whichpar], upper=upper[whichpar],method="Port")))
summary(Fit)             # not converged...

print(system.time(Fit<-modFit(p=parms(mc2)[whichpar],
                  f=Cost,lower=lower[whichpar], upper=upper[whichpar],method="L-BFGS-B")))
summary(Fit)             # not converged...

print(system.time(Fit<-modFit(p=parms(mc2)[whichpar],
                  f=Cost,lower=lower[whichpar], upper=upper[whichpar],method="BFGS")))
summary(Fit)             # not converged...

FITmrq <- modFit(Cost,parms(mc2)[whichpar],
                  lower=lower[whichpar], upper=upper[whichpar],method="Marq")
summary(FITmrq)

# 2. nls.lm fits the model to the data, using the residuals
FitMrq<- modFit(p=parms(mc2)[whichpar],f=Cost,method="Marq")
FitMrq[]
summary(FitMrq)

# a dirty trick, in case covariance matrix is singular
diag(FitMrq$hessian) <- diag(FitMrq$hessian) + 1e-11
summary(FitMrq)

BestPar <- FitMrq$par
plot(Cost(Fit$par),legpos="bottomleft")
Cost(FitMrq$par)
Cost(Oderes$par)

#===========================
# 3. Identifiability
#===========================
# Function that returns the modeled values at datapoints,
# as a function of input parameters
idFun <- function (p) Cost(p)$residual$mod

# Sensitivity functions at "best" values
sF<-sensFun( func=idFun,parms=BestPar,map=NULL)
summary(sF)

pairs(sF)

collin(sF)
# small enough value for all parameters
collin(sF,parset=1:length(whichpar))

# or, equivalently - and quicker!
sF2 <- sensFun(func=Cost,parms=BestPar)
collin(sF2,1:length(whichpar))  # collinearity of all parameters...

#=====================
# 4. MCMC application
#=====================

CM <- Cost(FitMrq$par)

# The true variance of the observed data are used.
s2prior <- CM$var$SSR/9

# using the estimated covariance does not work...
sP <- summary(FitMrq)
Covar  <- sP$cov.scaled * (2.4^2)/6

# so, 20% of parameter value is used for initial jump length
MCMC <- modMCMC(f=Cost,p=FitMrq$par,jump=FitMrq$par*0.2,niter=5000,
                var0=s2prior,wvar0=1,updatecov=10,lower=rep(0, length(whichpar)),
                upper = c(cells=1e8, mcyst=1e9,p=500,mu=2,K=1e9,dM=1)[whichpar])

plot(MCMC,Full=TRUE)     # not quite converged
hist(MCMC,Full=TRUE,col="blue")
pairs(MCMC)
cor(MCMC$pars)


sR<-sensRange(parInput=MCMC$pars,func=Solver)
plot(summary(sR),what="cells",legpos="bottomright")
points(dat$time,dat$cells)

plot(summary(sR),what="mcyst",legpos=NULL)
points(dat$time,dat$mcyst)

