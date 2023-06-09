\name{modCost}
\alias{modCost}

\title{
  Calculates the Discrepancy of a Model Solution with Observations
}

\description{
  Given a solution of a model and observed data, estimates the
  residuals, and the variable and model costs (sum of squared residuals).
}
\usage{
modCost(model, obs, x = "time", y = NULL, err = NULL,
        weight = "none", scaleVar = FALSE, cost = NULL,  ...)
}

\arguments{
  \item{model }{model output, as generated by the integration routine or
    the steady-state solver, a matrix or a data.frame, with one column
    per dependent and independent variable.
  }
  \item{obs }{the observed data, either in long (database) format (name,
    x, y), a data.frame, or in wide (crosstable, or matrix) format - see
    details.
  }
  \item{x }{the name of the \emph{independent} variable; it should be a
    name occurring both in the \code{obs} and \code{model} data
    structures.
  }
  \item{y }{either \code{NULL}, the name of the column with the
    \emph{dependent} variable values,or an index to the dependent
    variable values; if \code{NULL} then the observations are assumed to
    be in crosstable (matrix) format, and the names of the independent
    variables are given by the column names of this matrix.
  }
  \item{err }{either \code{NULL}, or the name of the column with the
    \emph{error} estimates, used to weigh the residuals (see details);
    if \code{NULL}, then the residuals are not weighed.
  }
  \item{cost }{if not \code{NULL}, the output of a previous call to
    \code{modCost}; in this case, the new output will combine both.
  }
  \item{weight}{only if \code{err}=\code{NULL}: how to weigh the
    residuals, one of "none", "std", "mean", see details.
  }
  \item{scaleVar}{if \code{TRUE}, then the residuals of one observed
    variable are scaled respectively to the number of observations (see
    details).
  }
  \item{...}{additional arguments passed to \R-function \code{approx}.
  }
}

\value{
  a list of type \code{modCost} containing:
  
  \item{model }{one value, the model cost, which equals the sum of scaled
    variable costs (see details).
  }
  \item{minlogp }{one value, -log(model probablity), where it is assumed
    that the data are normally distributed, with standard deviation =
    \code{error}.
  }
  \item{var}{ the variable costs, a data.frame with, for each observed
    variable the following (see details):
    \itemize{
      \item name, the name of the observed variable.
      \item scale, the scale-factor used to weigh the variable cost,
        either 1 or 1/(number observations), defaults to 1.
      \item N, the number of data points per observed variable.
      \item SSR.unweighted, the sum of squared residuals per observed
        variable, unweighted.
      \item SSR, the sum of weighted squared residuals per observed
        variable(see details).
    }
  }
  \item{residuals }{the data residual, a data.frame with several columns:
    \itemize{
      \item name, the name of the observed variable.
      \item x, the value of the independent variable (if present).
      \item obs, the observed variable value.
      \item mod, the corresponding modeled value.
      \item weight, the factor used to weigh the residuals, 1/error,
        defaults to 1.
      \item res, the weighted residuals between model and observations
        (mod-obs)*weight.
      \item res.unweighted, the residuals between model and observations
        (mod-obs).
    }
  }
}

\details{
  This function compares model output with observed data.

  It computes
  \enumerate{
    \item the weighted \emph{residuals}, one for each data point.
    \item the \emph{variable costs}, i.e. the sum of squared weight
      residuals per variable.
    \item the \emph{model cost}, the scaled sum of variable costs .
  }
  There are three steps:

  \emph{1.} For any observed data point, i, the \emph{weighted
  residuals} are estimated as:
  \deqn{res_i=\frac{Mod_i-Obs_i}{error_i}}{res_i=(mod_i-obs_i)/err_i}
  with \eqn{weight_i = 1/err_i} and where \eqn{Mod_i}{Mod_i} and
  \eqn{Obs_i}{Obs_i} are the modeled, respectively observed value of
  data point i.

  The weights are equal to 1/error, where the latter can be inputted,
  one for each data point by specifying \code{err} as an extra column in
  the observed data.

  This can only be done when the data input is in long (database) format.

  When \code{err} is not inputted, then the weights are specified via argument
  \code{weight} which is either:
  \itemize{
    \item \code{"none"}, which sets the weight equal to 1 (the default)
    \item \code{"std"}, which sets the weights equal to the reciprocal of the standard deviation
      of the observed data (can only be used if there is more than 1 data point)
    \item \code{"mean"}, which uses 1/mean of the absolute value of the
      observed data (can only be used if not 0).
  }
  \emph{2.} Then for each observed variable, j, a \emph{variable cost} is
  estimated as the sum of squared weighted residuals for this variable:
  
  \deqn{Cost_{var_j}= \sum \limits_{i = 1}^{n_j} {res_i}^2}{Cost_varj=sum(for i=1,n_j) (res_i^2)}

  where \eqn{n_j}{n_j} is the number of observations for observed
  variable j.

  \emph{3.} Finally, the \emph{model Cost} is estimated as the scaled sum of
  variable costs:

  \deqn{ModCost=\sum \limits_{j = 1}^{n_v} {\frac{Cost_{var_j}}{scale_{var_j}}}}{sum(Cost_varj/scale_j)}

  and where \eqn{scale_{var_j}}{scale_j} allows to scale the variable
  costs relative to the number of observations.  This is set by
  specifying argument \code{scaleVar}. If \code{TRUE}, then the variable
  costs are rescaled. The default is NOT to rescale
  (i.e. \eqn{scale_{var_j}}{scale_j}=1).

  The models typically consist of (a system of) differential equations, which
  are either solved by:
  \itemize{
    \item integration routines, e.g. the routines from package \code{deSolve},
    \item steady-state estimators, as from package \code{rootSolve}.
  }
  The data can be presented in two formats:
  \itemize{
    \item \emph{data table (long) format}; this is a two to four column
      data.frame that contains the \code{name} of the observed variable (always
      the FIRST column), the (optional) \code{value of the independent variable}
      (default column name = "time"), the \code{value of the observation} and
      the (optional) \code{value of the error}.
      For data presented in this format, the names of the column(s) with the
      independent variable (\code{x}) and the name of the column that has the
      value of the dependent variable \code{y} must be passed to function
      \code{modCost}.


    \item \emph{crosstable (wide) format}; this is a matrix, where each
      column denotes one dependent (or independent) variable; the column name
      is the name of the observed variable.
      When using this format, only the name of the column that contains the
      dependent variable must be specified (\code{x}).
  }
  As an example of both formats consider the data, called \code{Dat} consisting
  of two observed variables, called "Obs1" and "Obs2", both containing two
  observations, at time 1 and 2:
  \tabular{llll}{
  name    \tab time   \tab   val \tab err \cr
  Obs1    \tab 1      \tab   50  \tab 5  \cr
  Obs1    \tab 2      \tab  150  \tab 15  \cr
  Obs2    \tab 1      \tab  1    \tab 0.1  \cr
  Obs2    \tab 2      \tab  2    \tab 0.2  \cr
  }
  for the long format and

  \tabular{lll}{
   time   \tab   Obs1 \tab Obs2 \cr
   1      \tab  50    \tab 1 \cr
   2      \tab  150   \tab 2 \cr
  }
  for the crosstab format. Note, that in the latter case it is not possible to
  provide separate errors per data point.

  By calling modCost several consecutive times (using the \code{cost} argument),
  it is possible to combine both types of data files.
}

\note{
  In the future, it should be possible to have more than one independent
  variable present. This is not yet implemented, but it should allow e.g.
  to fit time series of spatially dependent variables.
}

\examples{

## =======================================================================
## Type 1 input:  name, time, value
## =======================================================================

## Create new data: two observed variables, "a", "b"
Data <- data.frame(name = c(rep("a", 4), rep("b", 4)),
                   time = c(1:4, 2:5), val = c(runif(4), 1:4))

## "a nonsense model"
Mod <- function (t, y, par) {
  da <- 0
  db <- 1
  return(list(c(da, db)))
}

out <- ode(y = c(a = 0.5, b = 0.5), times = 0:6, func = Mod, parms = NULL)

Data   # Show
out

## The cost function
modCost(model = out, obs = Data, y = "val")

## The cost function with a data error added
Dat2 <- cbind(Data, Err = Data$val*0.1)  # error = 10% of value
modCost(model = out, obs = Dat2, y = "val", err = "Err")


## =======================================================================
## Type 2 input:  Matrix format; column names = variable names
## =======================================================================

## logistic growth model
TT    <- seq(1, 100, 2.5)
N0    <- 0.1
r     <- 0.5
K     <- 100

## analytical solution
Ana <- cbind(time = TT, N = K/(1 + (K/N0 - 1) * exp(-r*TT)))

## numeric solution
logist <- function(t, x, parms) {
  with(as.list(parms), {
    dx <- r * x[1] * (1 - x[1]/K)
    list(dx)
  })
}

time  <- 0:100
parms <- c(r = r, K = K)
x     <- c(N = N0)

## Compare several numerical solutions
Euler <- ode(x, time, logist, parms, hini = 2, method = "euler")
Rk4   <- ode(x, time, logist, parms, hini = 2, method = "rk4")
Lsoda <- ode(x, time, logist, parms) # lsoda is default method
Ana2  <- cbind(time = time, N = K/(1 + (K/N0 - 1) * exp(-r * time)))

## the SSR and residuals with respect to the "data"
cEuler <- modCost(Euler, Ana)$model
cRk4   <- modCost(Rk4  , Ana)$model
cLsoda <- modCost(Lsoda, Ana)$model
cAna   <- modCost(Ana2 , Ana)$model
compare <- data.frame(method = c("euler", "rk4", "lsoda", "Ana"),
                      cost   = c(cEuler, cRk4, cLsoda, cAna))
## Plot Euler, RK and analytic solution
plot(Euler, Rk4, col = c("red", "blue"), obs = Ana,  
     main = "logistic growth", xlab = "time", ylab = "N")

legend("bottomright", c("exact", "euler", "rk4"), pch = c(1, NA, NA),
       col = c("black", "red", "blue"), lty = c(NA, 1, 2))
legend("right", ncol = 2, title = "SSR",
       legend = c(as.character(compare[,1]), 
                  format(compare[,2], digits = 2)))

compare

## =======================================================================
## Now suppose we do not know K and r and they are to be fitted...
## The "observations" are the analytical solution
## =======================================================================

## Run the model with initial guess: K = 10, r = 2
parms["K"] <- 10
parms["r"] <-  2

init <-  ode(x, time, logist, parms)

## FITTING algorithm uses modFit
## First define the objective function (model cost) to be minimised

## more general: using modFit
Cost <- function(P) {
  parms["K"] <- P[1]
  parms["r"] <- P[2]
  out <- ode(x, time, logist, parms)
  return(modCost(out, Ana))
}
(Fit<-modFit(p = c(K = 10, r = 2), f = Cost))

summary(Fit)

## run model with the optimized value:
parms[c("K", "r")] <- Fit$par
fitted <-  ode(x, time, logist, parms)

## show results, compared with "observations"
plot(init, fitted, col = c("green", "blue"), lwd = 2, lty = 1, 
     obs = Ana, obspar = list(col = "red", pch = 16, cex = 2), 
     main = "logistic growth", xlab = "time", ylab = "N")

legend("right", c("initial", "fitted"), col = c("green", "blue"), lwd = 2)

Cost(Fit$par)

}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>}

\references{
  Soetaert, K. and Petzoldt, T. 2010.  Inverse Modelling, Sensitivity and
  Monte Carlo Analysis in R Using Package FME.  Journal of Statistical
  Software 33(3) 1--28. \doi{10.18637/jss.v033.i03}
}
\keyword{utilities}

