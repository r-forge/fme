\name{sensRange}
\alias{sensRange}
\alias{summary.sensRange}
\alias{plot.sensRange}
\alias{plot.summary.sensRange}

\title{Global sensitivity analysis}
\description{Given a model consisting of differential equations, estimates
the global effect of certain (sensitivity) parameters on selected sensitivity variables.

This is done by drawing parameter values according to some predefined distribution, running
the model with each of these parameter combinations, and calculating summary statistics for
the selected output variables at each output interval.

This function thus produces 'envelopes' around the sensitivity variables.

}
\usage{
sensRange(func,parms,sensvar=NULL,dist="unif",
 parRange=NULL, parMean=NULL, parCovar=NULL, map=1,
 num=100,  ...)

\method{summary}{sensRange}(object,  \dots)
\method{plot}{summary.sensRange}(x, main=NULL, legpos="topleft", \dots)
\method{plot}{sensRange}(x, main=NULL, \dots)

 }
\arguments{
  \item{func }{An R-function that has as first argument \code{parms} and that returns a matrix with
  the values of the output variables (columns), and -optionally- a mapping variable (by default the first column)}
  \item{parms }{parameters passed to \code{func}; should be either a vector, or a list with named elements
    or a matrix of dimension (num, npar) with the values of the sensitivity parameters}
  \item{sensvar }{the output variables for which the sensitivity needs to be estimated. Either \code{NULL}, the default=all output variables,
  or a vector with output variable \code{names} (which should be present in the matrix returned by \code{func}),
  or a vector with \code{indices} to output variables as present in the output matrix
  (note that the column of this matrix with the mapping variable should not be selected) }
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
  \item{num }{the number of times the model has to be run. Set large enough. If \code{parms} is a matrix,
  then \code{num} parameters are selected randomly (from the rows of \code{parms}}
  \item{map }{The column number with the (independent) mapping variable in the output matrix returned by \code{func}. For
  dynamic models solved by integration, this will be the (first) column with \code{time}.
  Set to NULL if there is no mapping variable. Mapping variables cannot be selected for estimating
  sensitivity ranges; they are used for plotting.}
  \item{object }{an object of class \code{sensRange}}
  \item{x }{an object of class \code{sensRange}}
  \item{legpos }{position of the legend; set to NULL to avoid legend}
  \item{main }{a main title for the plot}
  \item{... }{additional arguments passed to function \code{func} or to the methods}
  }
\value{
  a data.frame of type \code{sensRange} containing the parameter set and
  the corresponding values of the sensitivity output variables.
  
  The list returned by \code{sensRange} has a method for the generic functions
   \code{\link{summary}},\code{\link{plot}} and \code{plot.summary}
 - see note.

  }

\details{

Models solved by integration (i.e. by using one of \code{'ode', 'ode.1D', 'ode.band', 'ode.2D'}),
have the output already in a form usable by \code{sensRange}.

The output for models solved by a steady-state solver (i.e. one of \code{steady, steady.1D,
steady.band, steady.2D}, needs to be rearranged - see examples.
}
\note{
The following \emph{methods} are included:

 \code{\link{summary}}, estimates summary statistics for the sensitivity variables,
 a data.frame with as many rows as there are mapping variables (or rows in the matrix returned by \code{func})
 and the following columns: \code{x}, the mapping value, \code{Mean}, the mean, \code{sd}, the standard deviation,
  \code{Min}, the minimal value, \code{Max}, the maximal value, \code{q25},  \code{q50}, \code{q75},
  the 25th, 50 and 75\% quantile
  
 \code{\link{plot}}, produces a "matplot" of the \code{sensRange} output, one plot for each sensitivity variable and
 with the mapping variable on the x-axis.

 \code{summary.plot}, produces a plot of the summary of the \code{sensRange} output, one plot for each sensitivity variable and
 with the ranges and mean +- standard deviation as coloured polygons.
  }
\examples{

#========================================
# Model from Soetaert and Herman, 2009
#========================================


pars <- list(gmax =0.5,eff = 0.5,
              ks =0.5, rB =0.01, dB =0.01)

solveBact <- function(pars)
{
 derivs <- function(t,state,pars)     # returns rate of change
  {
  with (as.list(c(state,pars)), {

  dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
  dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact

  return(list(c(dBact,dSub)))
                              })
  }

 state   <- c(Bact=0.1,Sub = 100)
 tout    <- seq(0,50,by=0.5)
 # ode solves the model by integration...
 return(as.data.frame(ode(y=state,times=tout,func=derivs,parms=pars)))
}

out <- solveBact(pars)

mf  <-par(mfrow=c(2,2))

plot(out$time,out$Bact,main="Bacteria",,
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)

# the sensitivity parameters
parRanges <- data.frame(min=c(0.4,0.4,0.),max=c(0.6,0.6,0.02))
rownames(parRanges)<- c("gmax","eff","rB")
parRanges

tout    <- 0:50
# sensitivity to rB; equally-spaced parameters ("grid")
SensR <- sensRange(func=solveBact,parms=pars,dist="grid",
                   sensvar="Bact",parRange=parRanges[3,],num=50)

Sens  <-summary(SensR)
plot(Sens,legpos="topleft",xlab="time, hour",ylab="molC/m3",
     main="Sensitivity to rB")

# sensitivity to all; latin hypercube
Sens2 <- summary(sensRange(func=solveBact,parms=pars,dist="latin",
                   sensvar=c("Bact","Sub"),parRange=parRanges,num=50))
plot(Sens2,xlab="time, hour",ylab="molC/m3",
     main="Sensitivity to gmax,eff,rB")

par(mfrow=mf )
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\keyword{utilities}

