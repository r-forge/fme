\name{sensFun}
\alias{sensFun}
\alias{summary.sensFun}
\alias{plot.sensFun}
\alias{pairs.sensFun}
\title{Local sensitivity analysis}
\description{Given a model consisting of differential equations, estimates
the local effect of certain (sensitivity) parameters on selected sensitivity variables
by calculating so-called sensitivity functions. This is a matrix whose (i,j)-th element contains
\deqn{ \frac{\partial y_i}{\partial \Theta _j}\cdot \frac{\Delta \Theta _j}{\Delta y_i}}{dy_i/dpar_j*parscale_j/varscale_i}
and where \eqn{y_i}{y_i} is an output variable (at a certain time instance), \eqn{\Theta_j}{par_j} is a parameter, and
 \eqn{\Delta y_i}{varscale_i} is the scaling of variable \eqn{y_i}{y_i},
 \eqn{\Delta \Theta_j}{parscale_j} is the scaling of parameter \eqn{\Theta_j}{par_j}.

}
\usage{
sensFun(func,parms,sensvar=NULL,senspar=names(parms),
 varscale=NULL, parscale=NULL,
 tiny=1e-8, map = 1, ...)

\method{summary}{sensFun}(object, vars=FALSE, \dots)
\method{pairs}{sensFun}(x,  \dots)
\method{plot}{sensFun}(x, main=NULL, legpos="topleft", \dots)

  }
\arguments{
  \item{func }{An R-function that has as first argument \code{parms} and that returns a matrix with
  the values of the output variables (columns), and -optionally- a mapping variable (by default the first column)}
  \item{parms }{parameters passed as the first argument to \code{solver}} or \code{func}
  \item{sensvar }{the output variables for which the sensitivity needs to be estimated. Either \code{NULL}, the default which will select all output variables,
  or a vector with output variable \code{names} (which should be present in model output matrix),
  or a vector with \code{indices} to output variables as present in the model output matrix
  }
  \item{senspar }{the parameters whose sensitivity needs to be estimated, the default=all parameters.
  Either a vector with parameter \code{names}, or a vector with \code{indices} to positions of parameters in \code{parms} }
  \item{varscale }{the scaling (weighing) factor for sensitivity variables, \code{NA}
                 indicates that the variable value is used. }
  \item{parscale }{the scaling (weighing) factor for sensitivity parameters, \code{NA}
                  indicates that the parameter value is used. }
  \item{tiny }{the perturbation, or numerical difference, factor, see details. }
  \item{map }{The column number with the (independent) mapping variable in the output matrix returned by \code{func}. For
  dynamic models solved by integration, this will be the (first) column with \code{time}.
  Set to NULL if there is no mapping variable. Mapping variables cannot be selected for estimating
  sensitivity functions; they are used for plotting.}
  \item{... }{additional arguments passed to \code{func} or to the methods}
  \item{object }{an object of class \code{sensFun}}
  \item{x }{an object of class \code{sensFun}}
  \item{vars }{if FALSE: summaries per parameter are returned; if TRUE: summaries per parameter and per variable are returned.}
  \item{legpos }{position of the legend; set to NULL to avoid legend}
  \item{main }{a main title for the plot}
  }
\value{
  a data.frame of class \code{sensFun} containing the sensitivity functions
  this is one row for each sensitivity variable at each
  independent (time or position) value and the following columns:
  \code{x}, the value of the independent (mapping) variable, usually time (solver= "ode.."), or distance (solver= "steady.1D")

  \code{var}, the name of the observed variable,

  \code{...}, a number of columns, one for each sensitivity parameter

  The data.frame returned by \code{sensFun} has methods for the generic functions
   \code{\link{summary}},  \code{\link{plot}},   \code{\link{pairs}}- see note.

}
\details{
There are essentially three different ways in which the model can be solved.

For each sensitivity parameter, the number of sensitivity functions estimated is:
length(sensvar) * length(mapping variable) and where \code{sensvar} are the selected sensitivity variables.


The sensitivity functions are estimated numerically. This means that
each parameter value \eqn{\Theta_j}{par_j} is perturbed as \eqn{\max{(tiny,\Theta_j *(1+tiny))}}{max(tiny,par_j)*(1+tiny)}
}
\note{
The following methods are provided:

\item \emph{summary}. Produces summary statistics of the sensitivity functions, a data.frame with:
 one row for each parameter and the following columns:

    \code{L1}: the L1-norm \eqn{\frac{1}{n} \cdot \sum{|S_{ij}|}}{sum(abs(Sij))/n}

    \code{L2}: the L2-norm \eqn{\frac{1}{n} \cdot \sqrt{\sum{S_{ij} \cdot S_{ij}}}}{sqrt(sum(Sij^2))/n}

    \code{Mean}: the mean of the sensitivity functions

    \code{Min}: the minimal value of the sensitivity functions

    \code{Max}: the maximal value of the sensitivity functions

  \item{var }{the summary of the variables sensitivity functions, a data.frame with the same columns as \code{model}
  and one row for each parameter + variable combination. This is only outputted if the variable names are effectively known}

\item \emph{plot } plots the sensitivity functions for each parameter

\item \emph{pairs } produces a pairs plot of the sensitivity results; overrides the default
  \code{gap = 0}, \code{upper.panel = NA}, and \code{diag.panel}

}
\examples{

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

plot(out$time,out$Bact,ylim=range(c(out$Bact,out$Sub)),
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)
lines(out$time,out$Sub,lty=2,lwd=2)
lines(out$time,out$Sub+out$Bact)

legend("topright",c("Bacteria","Glucose","TOC"),
       lty=c(1,2,1),lwd=c(2,2,1))

# sensitivity functions
SnsBact <- sensFun(func=solveBact,parms=pars,
                   sensvar="Bact",varscale=1)
head(SnsBact)
plot(SnsBact)
summary(SnsBact)

SF<- sensFun(func=solveBact,parms=pars,
             sensvar=c("Bact","Sub"),varscale=1)
head(SF)
tail(SF)

summary(SF,var=TRUE)

# Bivariate sensitivity
pairs(SnsBact)
mtext(outer=TRUE,side=3,line=-2,
      "Sensitivity functions",cex=1.5)

# pairwise correlation
cor(SnsBact[,-(1:2)])

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\references{
Soetaert, K and P.M.J. Herman, 2009. A practical guide to ecological modelling -
using R as a simulation platform. Springer, 390 pp.

Brun, R., Reichert, P., Kunsch, H.R., 2001.
Practical identificability analysis of large environmental simulation models.
Water Resour. Res. 37(4): 1015-1030
}
\keyword{utilities}

