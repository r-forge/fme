\name{sensRange}
\alias{sensRange}
\title{Global sensitivity analysis}
\description{Given a model consisting of differential equations, estimates
the global effect of certain (sensitivity) parameters on selected sensitivity variables.

This is done by drawing parameter values according to some predefined distribution, running
the model with each of these parameter combinations, and calculating summary statistics for
the selected output variables at each output interval.

This function thus produces 'envelopes' around the sensitivity variables.

}
\usage{
sensRange(parms,sensvar=NULL,dist="unif",
 parRange=NULL, parMean=NULL, parCovar=NULL, solver="ode",
 num=100, Full=FALSE, ...)
 }
\arguments{
  \item{parms }{parameters passed to the solver; should be either a vector, or a list with named elements
    or a matrix of dimension (num, npar) with the values of the sensitivity parameters}
  \item{sensvar }{the output variables for which the sensitivity needs to be estimated. Either \code{NULL}, the default=all output variables,
  or a vector with output variable \code{names} (which should be present in model output matrix),
  or a vector with \code{indices} to output variables as present in the model output matrix
  (note that the first column of this matrix is "times" and should not be selected) }
  \item{dist }{the distribution according to which parameter should be generated, one of
  "unif" (uniformly random samples), "norm", (normally distributed random samples),
  "latin" (latin hypercube distribution), "grid" (parameters arranged on a grid).
  The parameter distribution is specified by argument \code{parRange}, except for the
  normally distributed parameters, in which case the distribution is specified by the
  parameter means \code{parMean} and variance-covariance matrix, \code{parCovar}. Note that,
  if the distribution is \code{"norm"} and \code{parRange} is given, then a truncated distribution
  will be generated. (This is useful to prevent for instance that certain parameters become negative).
  Ignored if \code{parms} is a matrix}
  \item{parRange }{the range (min, max) of the sensitivity parameters, a matrix or (preferred) a data.frame with
  one row for each parameter, and two columns with the minimum (1st) and maximum (2nd) value.
  The rownames of \code{parRange} should be parameter names that are known in argument \code{parms}
  Ignored if \code{parms} is a matrix
  }
  \item{parMean }{only when \code{dist} is \code{"norm"}: the mean value of each parameter   }
  \item{parCovar }{only when \code{dist} is \code{"norm"}: the parameters variance-covariance matrix  }
  \item{solver }{the name of the solver routine, one of "ode", "ode.1D", "ode.2D", "steady", "steady.1D", "steady.2D"}
  \item{num }{the number of times the model has to be run. Set large enough. If \code{parms} is a matrix,
  then \code{num} parameters are selected randomly (from the rows of \code{parms}}
  \item{Full }{when TRUE full output is given, including all model runs and the parameter samples, when FALSE, only the summaries are returned}
  \item{... }{additional arguments passed to the integrator \code{ode}}
  }
\value{
  a list containing:
  \item{summ }{the summary of sensitivity variables, a data.frame with as many rows as there are
  sensitivity output
  variables (see details) and the following columns: \code{Mean}, the mean, \code{sd}, the standard deviation,
  \code{Min}{, the minimal value, \code{Max}, the maximal value, \code{q25},  \code{q50}, \code{q75},
  the 25th, 50 and 75\% quantile}
  }
  \item{sens }{Only when \code{Full} is TRUE: the full output, consisting of the parameter set and
  the corresponding values of the sensitivity output variables.}
}
\details{
When the model consists of \emph{differential equations} that are solved either by integration (package \code{deSolve} or
by steady-state calculation (package \code{rootSolve}), then the following (named) input arguments should be included:

  \item{func }{an \R-function that computes the values of the derivatives in the ODE system (the model definition) at time t.
       \code{func} must be defined as:
      \code{yprime = func(t, y, parms,...)}.  \code{t} is the current time point
      in the integration, \code{y} is the current estimate of the variables
      in the ODE system.  See details of integration routines in package \code{deSolve}}
  \item{y }{the initial (state) values for the ODE system, a vector. }
  \item{times }{if integrated: time sequence for which output is wanted; the first value of \code{times} must be the initial time}
  \item{time }{if steady-state: time at which steady-state should be estimated}


When the model is solved by integration (i.e. \code{solver} is one of \code{'ode', 'ode.1D', 'ode.band', 'ode.2D'}),
then for each sensitivity variable, as specified through \code{sensvar}, several output values are generated, one
for each element in \code{times}.

When the steady-state of the model is solved (i.e. \code{solver} is one of \code{steady, steady.1D,
steady.band, steady.2D}, then only one output value is generated.
}
\examples{

#========================================
# Model from Soetaert and Herman, 2009
#========================================


pars <- list(gmax =0.5,eff = 0.5,
              ks =0.5, rB =0.01, dB =0.01)

model <- function(t,state,pars)
{
with (as.list(c(state,pars)), {

dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact

return(list(c(dBact,dSub)))
                              })
}

Bini=0.1
Sini=100

state   <- c(Bact=Bini,Sub = Sini)
tout    <- seq(0,50,by=0.5)
out     <- as.data.frame(ode(y=state,times=tout,func=model,parms=pars))
mf<-par(mfrow=c(2,2))

plot(out$time,out$Bact,main="Bacteria",,
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)

# the sensitivity parameters
parRanges <- data.frame(min=c(0.4,0.4,0.),max=c(0.6,0.6,0.02))
rownames(parRanges)<- c("gmax","eff","rB")
parRanges

tout    <- 0:50
# sensitivity to rB; equally-spaced parameters ("grid")
Sens <- sensRange(func=model,y=state,times=tout,parms=pars,dist="grid",
                   sensvar="Bact",parRange=parRanges[3,],num=50)$summ
yrange<-range(cbind(Sens$Min,Sens$Max))
plot(tout,ylim=yrange,xlab="time, hour",ylab="molC/m3",type="n",
     main="Sensitivity rB")

polygon(c(tout,rev(tout)),c(Sens$Min,rev(Sens$Max)),
        col=grey(0.9),border=NA)
polygon(c(tout,rev(tout)),c(Sens$Mean-Sens$Sd,
        rev(Sens$Mean+Sens$Sd)),col=grey(0.8),border=NA)
lines(tout,Sens$Mean,lwd=2)
legend("topleft",fill=c(grey(0.9),grey(0.8)),
       legend=c("Min-Max","Mean+-sd"),bty="n")
legend("left",lty=1,lwd=2,legend="Mean",bty="n")

# sensitivity to all; latin hypercube
Sens2 <- sensRange(func=model,y=state,times=tout,parms=pars,dist="latin",
                   sensvar="Bact",parRange=parRanges,num=50)$summ
yrange<-range(cbind(Sens2$Min,Sens2$Max))
plot(tout,ylim=yrange,xlab="time, hour",ylab="molC/m3",type="n",
     main="Sensitivity gmax,eff,rB")

polygon(c(tout,rev(tout)),c(Sens2$Min,rev(Sens2$Max)),
        col=grey(0.9),border=NA)
polygon(c(tout,rev(tout)),c(Sens2$Mean-Sens2$Sd,
          rev(Sens2$Mean+Sens2$Sd)),col=grey(0.8),border=NA)
lines(tout,Sens2$Mean,lwd=2)

par(mfrow=mf )
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\keyword{utilities}

