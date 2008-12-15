\name{collin}
\alias{collin}
\title{Given the sensitivity functions of a model, estimates the collinearity of parameter sets}
\description{Based on the sensitivity functions of model variables to selection of parameters,
calculates the identifiability of sets of parameter.

The sensitivity functions are a matrix whose(i,j)-th element contains
\deqn{ \frac{\partial y_i}{\partial \Theta _j}\cdot \frac{\Delta (\Theta _j)}{\Delta y_i}}{dy_i/dpar_j*parscale_j/varscale_i}
and where \eqn{y_i}{y_i} is an output variable (at a certain time instance),
 \eqn{\Delta y_i}{varscale_i} is the scaling of variable \eqn{y_i}{y_i},
 \eqn{\Delta \Theta_j}{parscale_j} is the scaling of parameter \eqn{\Theta_j}{par_j}.

Function \code{collin} estimates the collinearity, or identifiability of all parameter sets.

In general a collinearity value less than about 20 is "identifiable"
}
\usage{collin(sensfun)}
\arguments{
  \item{sensfun }{Model sensitivity functions as estimated by \code{SensFun}}
  }
\value{
  a matrix with one row for each parameter combination (parameters as in \code{sensfun})

  Each row contains:

  for each parameter whether it is present (1) or absent (0) in the set,

  the number of parameters in the set and

  the collinearity value.
}

\details{
The collinearity is a measure of approximate linear dependence
between sets of parameters. The higher its value, the more the paramters are related.

It is possible to use the parameter identifiability technique for selecting
parameter sets that can be fine-tuned based on a data set.

Thus it is a powerful technique to make model calibration routines more robust,
because calibration routines often fail when parameters are strongly related.
In general, when the collinearity index exceeds 20, the linear dependence is assumed to be critical
(i.e. it will not be possible or easy to estimate all the parameters in the combination together).

The procedure is explained in Omlin et al. (2001).

1. First the parameter identifiability technique, is used to test how far a dataset
 can be used for estimating certain (combinaties of) parameters.
 After selection of these 'identifiable parameters' they are fine-tuned by calibration.
 
2. As the sensitivity analysis is a local analysis (i.e. its outcome depends on the current
 values of the model parameters) and the calibration routine is used to estimate the best
 values of the parameters, this is an iterative procedure.
 This means that identifiable parameters are determined, fitted to the data,
 then determined again, fitted, etcetera until convergenc is reached.
}
\examples{

#========================================
# Bacterial model from Soetaert and Herman, 2009
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

SnsBact <- sensFun(func=model,y=state,times=tout,parms=pars,
                   sensvar="Bact",varscale=1)$fun
matplot(tout,SnsBact[,-(1:2)],type="l")

#========================================
# Collinearity
#========================================
Coll <- collin(SnsBact[,-(1:2)])

# The larger the collinearity, the less identifiable the data set
format(data.frame(Coll),digits=2,scientific=FALSE)

nc <- ncol(Coll)
plot(Coll[,nc-1],Coll[,nc],log="y",main="Collinearity",
     xlab="Nr parameters",ylab="coll")

# 20 = magical number above which there are identifiability problems
abline(h=20,col="red")
Coll [Coll[,"collinearity"]<20 & Coll[,"N"]==4,]

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\references{
Soetaert, K and P.M.J. Herman, 2009. A practical guide to ecological modelling:
using R as a simulation platform. Springer, 390 pp.

Brun, R., Reichert, P., Kunsch, H.R., 2001.
Practical identifiability analysis of large environmental simulation models. Water Resour. Res. 37(4): 1015 1030

Omlin, M., Brun, R. and Reichert, P., 2001. Biogeochemical model of Lake Zurich:
sensitivity, identifiability and uncertainty analysis. Ecol. Modell. 141: 105 123.

}
\keyword{utilities}

