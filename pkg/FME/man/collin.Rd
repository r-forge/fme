\name{collin}
\alias{collin}
\alias{plot.collin}
\alias{print.collin}
\title{Given the sensitivity functions of a model, estimates the collinearity of
parameter sets}
\description{Based on the sensitivity functions of model variables to a selection
of parameters, calculates the "identifiability" of sets of parameter.

The sensitivity functions are a matrix whose(i,j)-th element contains
\deqn{ \frac{\partial y_i}{\partial \Theta _j}\cdot \frac{\Delta \Theta _j}{\Delta y_i}}{dy_i/dpar_j*parscale_j/varscale_i}
and where \eqn{y_i}{y_i} is an output variable, at a certain (time) instance, i,
 \eqn{\Delta y_i}{varscale_i} is the scaling of variable \eqn{y_i}{y_i},
 \eqn{\Delta \Theta_j}{parscale_j} is the scaling of parameter \eqn{\Theta_j}{par_j}.

Function \code{collin} estimates the collinearity, or identifiability of all parameter sets or
of one parameter set.

In general a collinearity value less than about 20 is "identifiable"
}
\usage{collin(sensfun,parset=NULL)

\method{print}{collin}(x, \ldots)
\method{plot}{collin}(x, \ldots)

}
\arguments{
  \item{sensfun }{Model sensitivity functions as estimated by \code{SensFun}}
  \item{parset }{One selected parameter combination, a vector with indices to
  the parameters }
  \item{x }{an object of class \code{collin}}
  \item{... }{additional arguments passed to the methods}
  }
\value{
  a data.frame of clas \code{collin} with one row for each parameter combination
  (parameters as in \code{sensfun})

  Each row contains:
  \item{... }{for each parameter whether it is present (1) or absent (0) in the set}
  \item{N}{ the number of parameters in the set}
  \item{collinearity}{ the collinearity value}

  The data.frame returned by \code{collin} has methods for the generic functions
   \code{\link{print}} and \code{\link{plot}}
}

\details{
The collinearity is a measure of approximate linear dependence
between sets of parameters. The higher its value, the more the parameters are related.
With "related" is meant that several paraemter combinations may produce similar
values of the output variables.
}
\note{
It is possible to use \code{collin} for selecting
parameter sets that can be fine-tuned based on a data set.
Thus it is a powerful technique to make model calibration routines more robust,
because calibration routines often fail when parameters are strongly related.
In general, when the collinearity index exceeds 20, the linear dependence is
assumed to be critical (i.e. it will not be possible or easy to estimate all the
parameters in the combination together).

The procedure is explained in Omlin et al. (2001).

1. First the function \code{collin} is used to test how far a dataset
 can be used for estimating certain (combinations of) parameters.
 After selection of an 'identifiable parameter set' (which has a low "collinearity")
 they are fine-tuned by calibration.
 
2. As the sensitivity analysis is a \emph{local} analysis (i.e. its outcome depends on the current
 values of the model parameters) and the calibration routine is used to estimate the best
 values of the parameters, this is an iterative procedure.
 This means that identifiable parameters are determined, fitted to the data,
 then a newly identifiable parameter set is determined, fitted, etcetera until convergenc is reached.

See the paper by Omlin et al. (2001) for more information.
}
\examples{

#========================================
# Bacterial model from Soetaert and Herman, 2009
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
SnsBact <- sensFun(func=solveBact,parms=pars,
                   sensvar="Bact",varscale=1)

#========================================
# Collinearity
#========================================
Coll <- collin(SnsBact)

# The larger the collinearity, the less identifiable the data set
Coll

plot(Coll,log="y")

# 20 = magical number above which there are identifiability problems
abline(h=20,col="red")

Coll [Coll[,"collinearity"]<20 & Coll[,"N"]==4,]

# collinearity of one parameter set
collin(SnsBact,c(1,3,5))
collin(SnsBact,1:5)

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\references{
Brun, R., Reichert, P., Kunsch, H.R., 2001.
Practical identifiability analysis of large environmental simulation models. Water Resour. Res. 37(4): 1015 1030

Omlin, M., Brun, R. and Reichert, P., 2001. Biogeochemical model of Lake Zurich:
sensitivity, identifiability and uncertainty analysis. Ecol. Modell. 141: 105 123.

}
\keyword{utilities}

