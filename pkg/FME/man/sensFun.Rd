\name{sensFun}
\alias{sensFun}
\title{Local sensitivity analysis}
\description{Given a model consisting of differential equations, estimates
the local effect of certain (sensitivity) parameters on selected sensitivity variables
by calculating so-called sensitivity functions. This is a matrix whose (i,j)-th element contains
\deqn{ \frac{\partial y_i}{\partial \Theta _j}\cdot \frac{\Delta (\Theta _j)}{\Delta y_i}}{dy_i/dpar_j*parscale_j/varscale_i}
and where \eqn{y_i}{y_i} is an output variable (at a certain time instance), \eqn{\Theta_j}{par_j} is a parameter, and
 \eqn{\Delta y_i}{varscale_i} is the scaling of variable \eqn{y_i}{y_i},
 \eqn{\Delta \Theta_j}{parscale_j} is the scaling of parameter \eqn{\Theta_j}{par_j}.

The function also returns several summaries of these sensitivity functions
}
\arguments{
  \item{parms }{parameters passed as the first argument to \code{func} or \code{solver}}
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
  \item{solver }{either the name of the solver routine, one of \code{"ode", "ode.1D", "ode.2D", "steady", "steady.1D", "steady.2D"}, or
                 a function that takes as first argument the parameter values \code{parms} - see details}
  \item{... }{additional arguments passed to the solver}
  }
\usage{
sensFun(parms,sensvar=NULL,senspar=names(parms),
 varscale=NULL, parscale=NULL,
 tiny=1e-8, solver="ode", ...)
  }
\value{
  a list containing:
  
  \item{model }{the summary of the model sensitivity functions, a data.frame with one row for each parameter and
                 the following columns:

    \code{L1}: the L1-norm \eqn{\frac{1}{n} \cdot \sum{|S_{ij}|}}{sum(abs(Sij))/n}

    \code{L2}: the L2-norm \eqn{\frac{1}{n} \cdot \sum{\sqrt{S_{ij} \cdot S_{ij}}}}{sum(Sij^2)/n}

    \code{Mean}: the mean of the sensitivity functions

    \code{Min}: the minimal value of the sensitivity functions

    \code{Max}: the maximal value of the sensitivity functions
    }
  \item{var }{the summary of the variables sensitivity functions, a data.frame with the same columns as \code{model}
  and one row for each parameter + variable combination. This is only outputted if the variable names are effectively known}

  \item{fun }{a data.frame with the sensitivity functions, this is one row for each sensitivity variable at each
  independent (time or position) value.
   and the following columns:
   
                  * \code{x}, the value of the independent variable, usually time (solver= "ode.."), or distance (solver= "steady.1D")

                  * \code{var}, the name of the observed variable,

                  * \code{...}, a number of columns, one for each sensitivity parameter
                  }
}
\details{
If the model is solved by integration (e.g. using a function from package \code{deSolve}), then
\code{solver} should be one of (\code{"ode", "ode.1D", "ode.2D"}).
In this case, for each sensitivity parameter, the number of sensitivity functions estimated is:
length(sensvar) * length(times) and where \code{sensvar} are the selected sensitivity variables and
\code{times} is the number of time instances at which the model has been solved.

If the model is solved to steady-state (using a function from package \code{rootSolve}),
then \code{solver} should be one of (\code{"steady", "steady.1D", "steady.2D"}).
In this case the number of sensitivity functions estimated equals the length of \code{sensvar}

If the model is solved otherwise, then
\code{solver} should be a function that takes as first argument the parameter values \code{parms}, and that
returns the ouput, a matrix with column names equal to model variables. Alternatively, the output can be a vector (with named values).
If the function returns a vector or a matrix with one row, then the number of sensitivity functions estimated equals the length of \code{sensvar},
else the number of sensitivity functions will equal number of rows * length of \code{sensvar}.

The sensitivity functions are estimated numerically. This means that
each parameter value \eqn{\Theta_j}{par_j} is perturbed as \eqn{\max{(tiny,\Theta_j *(1+tiny))}}{max(tiny,par_j)*(1+tiny)}

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

tout    <- seq(0,50,by=0.5)
state   <- c(Bact=Bini,Sub = Sini)

out     <- as.data.frame(ode(state,tout,model,pars))

plot(out$time,out$Bact,ylim=range(c(out$Bact,out$Sub)),
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)
lines(out$time,out$Sub,lty=2,lwd=2)
lines(out$time,out$Sub+out$Bact)

legend("topright",c("Bacteria","Glucose","TOC"),
       lty=c(1,2,1),lwd=c(2,2,1))

SnsBact <- sensFun(func=model,y=state,times=tout,parms=pars,
                   sensvar="Bact",varscale=1)$fun
matplot(tout,SnsBact[,-(1:2)],type="l",ylab="Sfun",
        main="sensitivity functions")
legend("topleft",names(SnsBact[-(1:2)]),col=1:5,lty=1:5)

sensFun(func=model,y=state,times=tout,parms=pars,
        sensvar="Bact",varscale=1)$model

SF<- sensFun(func=model,y=state,times=tout,parms=pars,
             sensvar=c("Bact","Sub"),varscale=1)
SF$model
SF$var
# Bivariate sensitivity

panel.cor <- function(x, y)
             text(x=mean(range(x)),y=mean(range(y)),
             labels=format(cor(x,y),digits=2))

pairs(SnsBact[,-(1:2)],upper.panel=panel.cor)
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

