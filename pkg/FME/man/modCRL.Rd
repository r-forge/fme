\name{modCRL}
\alias{modCRL}
\alias{summary.modCRL}
\alias{plot.modCRL}
\title{
  Monte Carlo analysis
}
\description{
  Given a model consisting of differential equations, estimates
  the global effect of certain (sensitivity) parameters on selected
  sensitivity variables.

  This is done by drawing parameter values according to some predefined
  distribution, running the model with each of these parameter combinations,
  and calculating the values of the selected output variables at each output
  interval.

  This function is useful for "what-if" scenarios.

  If the output variables consist of a time-series or spatially dependent, use
  sensRange instead.

}
\usage{
modCRL(func,parms=NULL,sensvar=NULL,dist="unif",
  parInput=NULL, parRange=NULL, parMean=NULL, parCovar=NULL,
  num=100,  ...)

\method{summary}{modCRL}(object,  \dots)

\method{plot}{modCRL}(x, what=NULL,xlab=NULL, ylab=NULL, \dots)
}
\arguments{
  \item{func }{an R-function that has as first argument \code{parms} and that
    returns a vector with variables whose sensitivity should be estimated.
  }
  \item{parms }{parameters passed to \code{func}; should be either a vector,
    or a list with named elements.
    If \code{NULL}, then the first element of \code{parInput} is taken.
  }
  \item{sensvar }{the output variables for which the sensitivity needs to be
    estimated. Either \code{NULL}, the default=all output variables,
    or a vector with output variable \code{names} (which should be present in
    the vector returned by \code{func}),
    or a vector with \code{indices} to output variables as present in the
    output vector returned by \code{func}.
  }
  \item{dist }{the distribution according to which the parameters should be
    generated, one of \code{"unif"} (uniformly random samples),
    \code{"norm"}, (normally distributed random samples),
    \code{"latin"} (latin hypercube distribution),
    \code{"grid"} (parameters arranged on a grid).
    The input parameters for the distribution are specified by \code{parRange}
    (min,max), except for the normally distributed parameters, in which case
    the distribution is specified by the parameter means \code{parMean} and
    the variance-covariance matrix, \code{parCovar}. Note that, if the
    distribution is \code{"norm"} and \code{parRange} is given, then a
    truncated distribution will be generated.
    (This is useful to prevent for instance that certain parameters become
    negative). Ignored if \code{parInput} is specified.
  }
  \item{parRange }{the range (min, max) of the sensitivity parameters, a matrix
    or (preferred) a data.frame with one row for each parameter, and two
    columns with the minimum (1st) and maximum (2nd) value.
    The rownames of \code{parRange} should be parameter names that are known
    in argument \code{parms}. Ignored if \code{parInput} is specified.
  }
  \item{parInput }{a matrix with dimension (*, npar) with the values of the
    sensitivity parameters.
  }
  \item{parMean }{only when \code{dist} is \code{"norm"}: the mean value of
    each parameter. Ignored if \code{parInput} is specified.
  }
  \item{parCovar }{only when \code{dist} is \code{"norm"}: the parameter's
    variance-covariance matrix.
  }
  \item{num }{the number of times the model has to be run. Set large enough.
    If \code{parInput} is specified, then \code{num} parameters are selected
    randomly (from the rows of \code{parInput}.
  }
  \item{object }{an object of class \code{modCRL}.
  }
  \item{x }{an object of class \code{modCRL}.
  }
  \item{what }{the name or the index to the variables that should be plotted.
  }
  \item{xlab }{X axis label to the plot.
  }
  \item{ylab }{Y axis label to the plot.
  }
  \item{... }{additional arguments passed to function \code{func} or to
    the methods.
  }
}

\value{
  a data.frame of type \code{modCRL} containing the parameter(s) and
  the corresponding values of the sensitivity output variables.
  
  The list returned by \code{modCRL} has a method for the generic functions
  \code{\link{summary}} and \code{\link{plot}}
 - see note.
  }

\note{
The following \emph{methods} are included:
  \itemize{

    \item \link{summary}, estimates summary statistics for the sensitivity
      variables, a table with as many rows as there are variables (or elements
      in the vector returned by \code{func}) and the following columns:
      \code{x}, the mapping value, \code{Mean}, the mean, \code{sd}, the
      standard deviation, \code{Min}, the minimal value, \code{Max}, the
      maximal value, \code{q25},  \code{q50}, \code{q75}, the 25th, 50 and
      75\% quantile
  
    \item \link{plot}, produces a plot of the \code{modCRL} output, one plot
      for each sensitivity variable and with the parameter value on the x-axis.
      This only works when there is only one parameter!

  }
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

plot(out$time,out$Bact,main="Bacteria",,
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)

# Function that returns the last value of the simulation
SF <- function (p)
{
  pars["eff"]<-p
  out <- solveBact(pars)
  return(out[nrow(out),2:3])
}

CRL <- modCRL(func=SF,parRange=c(0.2,0.8))

plot(CRL)  # plots both variables
plot(CRL, what="Bact") #selects one

}
\author{
Karline Soetaert <k.soetaert@nioo.knaw.nl>.
}
\keyword{utilities}

