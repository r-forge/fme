\name{modFit}
\alias{modFit}
\alias{summary.modFit}
\alias{deviance.modFit}
\alias{coef.modFit}
\alias{residuals.modFit}
\alias{df.residual.modFit}
\alias{plot.modFit}
\alias{print.summary.modFit}
\title{
  Constrained Fitting of a Model to Data
}

\description{
  Fitting a model to data, with lower and/or upper bounds
}

\usage{modFit(f, p, ..., lower=-Inf, upper=Inf,
  method=c("Marq","Port","Newton",
           "Nelder-Mead","BFGS","CG","L-BFGS-B","SANN",
           "Pseudo"),
  control=list(),hessian=TRUE)

\method{summary}{modFit}(object, cov=TRUE,\dots)

\method{deviance}{modFit}(object,  \dots)

\method{coef}{modFit}(object, \dots)

\method{residuals}{modFit}(object, \dots)

\method{df.residual}{modFit}(object, \dots)

\method{plot}{modFit}(x, \dots)

\method{print}{summary.modFit}(x, digits = max(3, getOption("digits") - 3),
               \dots)
}
\arguments{
  \item{f }{a function to be minimized, with first argument the vector of
    parameters over which minimization is to take place. It should return
    either a vector of \emph{residuals} (of model versus data) or an element
    of class \emph{modCost} (as returned by a call to \code{\link{modCost}}.
  }
  \item{p }{initial values for the parameters to be optimized over.
  }
  \item{...}{additional arguments passed to function \code{f} (modFit) or passed
    to the methods.
  }
  \item{lower }{lower bounds on the parameters; if unbounded set equal to
    \code{-Inf}.
  }
  \item{upper }{upper bounds on the parameters; if unbounded set equal to
    \code{Inf}.
  }
  \item{method }{The method to be used, one of "Marq","Port","Newton",
    "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Pseudo" - see details.
  }
  \item{hessian }{TRUE if hessian is to be estimated. Note that, if set to
    \code{FALSE}, then a summary cannot be estimated.
  }
  \item{control }{additional control arguments passed to the optimisation
    routine - see details of \code{\link[minpack.lm]{nls.lm}} ("Marq"),
    \code{\link{nlminb}} ("Port"),
    \code{\link{optim}} ("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
    \code{\link{nlm}} ("Newton") or \code{\link{pseudoOptim}}("Pseudo").
  }
  \item{object }{an object of class \code{modFit}.
  }
  \item{x }{an object of class \code{modFit}.
  }
  \item{digits }{number of significant digits in printout.
  }
  \item{cov }{when \code{TRUE} also calculates the parameter covariances.
  }
}

\value{
  a list of class \emph{modFit} containing the results as returned from the
  called optimisation routines.
  
  This includes the following:
  
  \item{par }{the best set of parameters found.
  }
  \item{ssr }{the sum of squared residuals, evaluated at the best set of
    parameters.
  }
  \item{hessian }{A symmetric matrix giving an estimate of the Hessian
    at the solution found. - see note.
  }
  \item{residuals }{the result of the last \code{f} evaluation; that is, the
    residuals.
  }
  \item{ms }{the mean squared residuals, i.e. \code{ssr/length(residuals)}.
  }
  \item{var_ms }{the weighted and scaled variable mean squared residuals,
    one value per observed variable; only when \code{f}
    returns an element of class \emph{modCost} ; \code{NA} otherwise.
  }
  \item{var_ms_unscaled }{the weighted, but not scaled variable mean
    squared residuals
  }
  \item{var_ms_unweighted }{the raw variable mean squared residuals,
    unscaled and unweighted.
  }
  \item{... }{any other arguments returned by the called optimisation routine.
  }

  Note: this means that some return arguments of the original functions are
  renamed.
  
  More specifically, "objective" and "counts" from routine  \code{\link{nlminb}}
  (method=Port) are renamed; "value" and "counts"; "niter" and "minimum"
   from routine  \code{\link[minpack.lm]{nls.lm}} (method=Marq)
   are renamed; "counts" and "value"; "minimum" and "estimate" from routine
   \code{\link{nlm}} (method=Newton)
   are renamed;

  The list returned by \code{modFit} has a method for the generic functions
  \code{\link{summary}},\code{\link{deviance}},\code{\link{coef}},
  \code{\link{residuals}},\code{\link{df.residual}} and \code{print.summary}
  - see note.
}

\details{

  Note that arguments after ... must be matched exactly.

  The method to be used is specified by argument \code{method} which can be

  One of the methods from function \code{\link{optim}}:

  \itemize{
    \item Nelder-Mead, the default from \code{optim}
    \item BFGS, a quasi-Newton method
    \item CG, a conjugate-gradient method
    \item L-BFGS-B, constrained quasi-Newton method
    \item SANN, method of simulated annealing
  }
  Or one of the following:
  \itemize{
    \item Marq, the Levenberg-Marquardt algorithm
     (\code{\link[minpack.lm]{nls.lm}} from package \code{minpack}) -
     the default. Note that this method is the only least squares method.
    \item Newton, a Newton-type algorithm (see \code{\link{nlm}})
    \item Port, the Port algorithm (see \code{\link{nlminb}})
    \item Pseudo, the pseudorandom-search algorithm
      (see \code{\link{pseudoOptim}}).
  }

  For difficult problems it may be efficient to perform some iterations with
  \code{Pseudo}, which will bring the algorithm near the vicinity of a
  (the) minimum, after which the default algorithm (\code{Marq}) is used to
  locate the exact minimum.

  The implementation for the routines from \code{\link{optim}} differs from
  \code{\link{constrOptim}} which implements an adaptive barrier algorithm and
  which allows a more flexible implementation of linear constraints.

  In \code{modFit}, bounds on parameters are imposed by a transformation of the
  parameters to be fitted.

  In case \emph{both lower and upper bounds} are specified, this is achieved by
  by a tangents and arc tangents transformation.

  This is, parameter values, p', generated by the optimisation routine, and
  which are located in the range [-Inf,Inf] are transformed, before they are
  passed to \code{f} as:

   \deqn{p = (upper + lower)/2 + (upper - lower) \cdot \arctan(p')/\pi}{%
        p = 1/2 * (upper + lower) + (upper - lower) * arctan(p')/pi}.

  which maps them into the interval [lower,upper]

  Before the optimisation routine is called, the original parameter values, as
  given by argument \code{p} are mapped from [lower,upper] to [-Inf,Inf] by:

  \deqn{p' = \tan(\pi/2 \cdot (2 p - upper - lower) / (upper - lower)) }{%
       p' = tan(pi/2 * (2 * p - upper - lower) / (upper - lower))}


  In case \emph{only lower or upper bounds} are specified, this is achieved
  by a log transformation and a corresponding exponential back transformation.

  In case parameters are transformed, (all methods), or in case the
  \code{method} \code{Port}, \code{Pseudo} or \code{Marq} is selected,
  the \emph{hessian} is approximated as \eqn{2*J^T*J}, where J is the Jacobian,
  estimated by finite differences.

  This ignores the second derivative terms, but this is reasonable if the method
  has truly converged to the minimum.

  Note that finite differences are not extremely precise.

  In case the Levenberg-Marquard method (\code{Marq}) is used, and parameters
  are not transformed, 0.5 times the Hessian of the least squares problem is
  returned by \code{nls.lm}, the original marquardt algorithm. To make it
  compatible, this value is multiplied with 2 and the TRUE hessian is thus
  returned by \code{modFit}.
}

\note{
  The \code{\link{summary}} \emph{method} is based on an estimate of the
  parameter covariance matrix.
  In computing the covariance matrix of the fitted parameters, the problem is
  treated as if it were a linear least squares problem, linearizing around
  the parameter values that minimize \eqn{Chi^2}.

  The covariance matrix is estimated as 1/(0.5*hessian)

  This computation relies on several things, i.e.

  \enumerate{
    \item the parameter values are located at the minimum (i.e. the fitting
      algorithm has converged).
    \item the observations yj are subject to independent errors whose
      variances are well estimated by 1 / (n - p) times the residual sum
      of squares (where n= number of data points, p = number of parameters).
    \item the model is not too nonlinear
  }
  This means that the estimated covariance (correlation) matrix and the
  confidence intervals derived from it may be worthless if the assumptions
  behind the covariance computation are invalid.

  If in doubt about the validity of the summary computations, use monte carlo
  fitting instead, or run a \code{\link{modMCMC}}.

  Other methods included are:
  \itemize{
    \item \link{deviance}, which returns the model deviance,
    \item \link{coef}, which extracts the values of the fitted parameters,
    \item \link{residuals},which extracts the model residuals,
    \item \link{df.residual} which returns the residual degrees of freedom
    \item print.summary, producing a nice printout of the summary.
  }
}
\seealso{
  \code{\link{constrOptim}} for constrained optimization.
}

\examples{

# logistic growth model
TT    <- seq(1,60,5)
N0    <- 0.1
r     <- 0.5
K     <- 100

## perturbed analytical solution
Data <- data.frame(time=TT,N=K/(1+(K/N0-1)*exp(-r*TT))*(1+rnorm(length(TT),sd=0.01)))

plot(TT, Data[,"N"], ylim=c(0, 120), pch=16, col="red",
     main = "logistic growth", xlab="time", ylab="N")


####################################
# Fitted with analytical solution  #
####################################
# initial "guess"
parms <- c(r = 2, K = 10, N0=5)

# analytical solution
model <- function(parms,time)
  with (as.list(parms),return(K/(1+(K/N0-1)*exp(-r*time))))

# run the model with initial guess and plot results
lines (TT, model(parms,TT) ,lwd=2,col="green")

# FITTING algorithm 1
ModelCost <- function(P)
{
 out <-model(P,TT)
 return(Data$N-out)  # residuals
}

(Fita <- modFit(f=ModelCost,p=parms))

times <- 0:60
lines(times,model(Fita$par,times),lwd=2, col="blue")
summary(Fita)

####################################
#  Fitted with numerical solution  #
####################################
## numeric solution
logist <- function(t, x, parms) {
  with(as.list(parms), {

    dx <- r * x[1] * (1 - x[1]/K)

    list(dx)
  })
}

# model cost,
ModelCost2 <- function(P)
{
 out <- ode(y=c(N=P[["N0"]]),func=logist,parms=P,times=c(0,TT))
 return(modCost(out,Data))  # object of class modCost
}

Fit <- modFit(f=ModelCost2,p=parms,lower=rep(0,3),upper=c(5,150,10))

out <- ode(y=c(N=Fit$par[["N0"]]),func=logist,parms=Fit$par,times=times)

lines(out,col="red",lty=2)
legend("right",c("data","original","fitted analytical","fitted numerical"),
       lty=c(NA,1,1,2),lwd=c(NA,2,2,1),col=c("red","green","blue","red"),
       pch=c(16,NA,NA,NA))
summary(Fit)
plot(residuals(Fit))
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>,

  Thomas Petzoldt <thomas.petzoldt@tu-dresden.de>
}
\references{
  Press, W. H.,  Teukolsky, S. A., Vetterling, W. T. and
  Flannery, B. P. (2007) Numerical Recipes in C. Cambridge
  University Press.

  \url{http://netlib.bell-labs.com/cm/cs/cstr/153.pdf}
}
\keyword{utilities}

