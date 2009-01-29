library("simecol")
require(FME)

fem2mu <- 1e+6 * 1e-15
mu2fem <- 1e+9

plot.mcyst <- function(ex, sc) {
  o <- out(sc)
  nf <- layout(matrix(c(1, 2), 2, 2, byrow=FALSE), c(3, 1), 1)
  opar <- par(mar=c(0,6.1,2,1), pch=16, cex=1, las=1)
  ylim_cells <- c(0, max(o$cells, na.omit(ex$cells)))
  ylim_mcyst  <- c(0, max(o$mcyst*fem2mu, na.omit(ex$mcyst)))
  plot(ex$time, ex$cells, axes=FALSE, col="forestgreen", ylab="", ylim=ylim_cells)
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
  lines(o$time, o$mcyst * fem2mu)
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
  	mu = 0.238,       # 1/d (Lab)
    K  = 31.5e6,      # Cells per ml
  	p  = 120,         # fg/cell
  	dM = 0.04,        # 1/d (Orr & Jones)
  	cells = 1.4e6,
  	mcyst = 8.6 * 1e7
  ),
  #inputs = matrix(data.frame(time=0:100, pH=rep(8.6, 101))),
  times = c(from=0, to=100, by=1),
  solver = "lsoda",
  initfunc = function(obj) {
    #print(parms(obj))
    init(obj)[c("cells", "mcyst")] <- parms(obj)[c("cells", "mcyst")]
    obj
  }
)

# Read the data...
dat <- read.table("expS108.txt", header=TRUE, sep="\t")
dat$mcyst[1] <- 0 ## replace NA

# Make the data compatible with model...
dat$mcyst <- dat$mcyst * mu2fem


# Simulate,
mc2 <- mcyst2
sc  <- sim(mc2)
plot.mcyst(dat, sc)

# Parameters to be fitted...
whichpar <- c("cells", "mcyst", "p", "mu", "K", "dM")
lower    <- c(cells = 1e5, mcyst = 0,   p = 0,   mu = 0.1, K = 10e6, dM = 0)
upper    <- c(cells = 5e6, mcyst = 5e8, p = 120, mu = 1,   K = 40e6, dM = 1)

#===========================
# 1. Fit, using fitOdeModel
#===========================

obstime   <- dat$time
times(sc) <- obstime

optdat <- data.frame(cells=dat$cells, mcyst=dat$mcyst )

# this works...
Oderes <- fitOdeModel(sc, whichpar=whichpar, obstime, optdat,  method="BFGS",
  debuglevel=0, lower = lower, upper = upper, control=list(trace=TRUE))

Oderes$par

# this does not work...
Oderes2 <- fitOdeModel(sc, whichpar=whichpar, obstime, optdat,  method="L-BFGS-B",
  debuglevel=0, lower = lower, upper = upper, control=list(trace=TRUE))

parms(sc)[whichpar] <- Oderes$par

times(sc) <-  c(from=0, to=100, by=1)
sc        <- sim(sc)
plot.mcyst(dat, sc)

#parms(sc)
#           mu            K            p           dM
# 2.386138e-01 3.154855e+07 1.200000e+02 3.987729e-02
#        cells        mcyst
# 1.424767e+06 8.627639e+07


#===========================
# 2. Fit using other algorithms ...
#===========================
# Solver function
Solver <- function (p)
{
 parms(sc)[whichpar]<- p
 out(sim(sc))
}

# Cost of model versus data
Cost <- function (p)
{
 Run <- Solver(p)
 return ( modCost(model=Run,obs=dat,weight="std"))
}

print(system.time(Fit<-modFit(p=parms(mc2)[whichpar],
                  f=Cost,lower=lower, upper=upper,method="Port")))
summary(Fit)             # not converged...

print(system.time(Fit<-modFit(p=parms(mc2)[whichpar],
                  f=Cost,lower=lower, upper=upper,method="L-BFGS-B")))
summary(Fit)             # not converged...

print(system.time(Fit<-modFit(p=parms(mc2)[whichpar],
                  f=Cost,lower=lower, upper=upper,method="BFGS")))
summary(Fit)             # not converged...

FITmrq <- modFit(Cost,parms(mc2)[whichpar],
                  lower=lower, upper=upper,method="Marq")
summary(FITmrq)

# 2. nls.lm fits the model to the data, using the residuals
FitMrq<- modFit(p=parms(mc2)[whichpar],f=Cost,method="Marq")
FitMrq[]
summary(FitMrq)

# a dirty trick
diag(FitMrq$hessian) <- diag(FitMrq$hessian) + 1e-11
summary(FitMrq)

BestPar <- FitMrq$par
plot(Cost(Fit$par))
Cost(FitMrq$par)
Cost(Oderes$par)

#===========================
# 3. Identifiability
#===========================
# Function that returns the modeled values at datpoints, as a function of input parameters
idFun <- function (p) Cost(p)$residual$mod

# Sensitivity functions at "best" values
sF<-sensFun( func=idFun,parms=BestPar,map=NULL)
summary(sF)

pairs(sF)

collin(sF)

# or, equivalently - and quicker!
sF2 <- sensFun(func=Cost,parms=BestPar)
collin(sF2,1:6)  # collinearity of all parameters...

#=====================
# 4. MCMC application
#=====================

CM <- Cost(FitMrq$par)

# The true variance of the observed data are used.
s2prior <- CM$var$SSR.unweighted/9

# using the estimated covariance does not work...
sP <- summary(FitMrq)
Covar  <- sP$cov.unscaled * (sum(s2prior))/2 * (2.4^2)/6

# so, 20% of parameter value is used for initial jump length
MCMC <- modMCMC(f=Cost,p=FitMrq$par,jump=FitMrq$par*0.2,niter=5000,
                var0=s2prior,wvar0=1,updatecov=10,lower=c(0,0,0,0,0,0),
                upper = c(1e8,1e9,500,2,1e9,1))

plot(MCMC,Full=TRUE)     # not quite converged
hist(MCMC,Full=TRUE,col="blue")
pairs(MCMC)
cor(MCMC$pars)


sR<-sensRange(parInput=MCMC$pars,func=Solver)
plot(summary(sR),what="cells",legpos="bottomright")
points(dat$time,dat$cells)

plot(summary(sR),what="mcyst",legpos=NULL)
points(dat$time,dat$mcyst)

