\name{FME-package}
\alias{FME-package}
\alias{FME}
\docType{package}
\title{
Flexible Modelling Environment - Utility functions for models consisting of
differential equations
}
\description{
R-package FME contains functions to run complex applications of models that
produce output as a function of input parameters. It is mainly intended to be
used with models consisting of ordinary differential equations (ODE),
partial differential equations (PDE) or differential algebraic equations (DAE).

It contains:
\itemize{
\item Functions to allow fitting of the model to data. Function \code{modCost}
estimates the residuals between model output and data, (where the data are
represented as a function of one or more independent variables) and estimates, a.o.
the sum of weighted squared residuals. Function \code{modFit} uses the output
of \code{modCost} to find the best-fit parameters. It provides a wrapper around
\R's built-in minimisation routines (\code{optim}, \code{nlm}, \code{nlminb})
and \code{nls.lm} from package \code{minpack.lm}.
FME also includes an implementation of the pseudo-random search algorithm
(function \code{pseudoOptim}).

\item Function \code{sensFun} estimates the sensitivity functions of selected
output variables as a function of model parameters. This is the basis of
uni-variate, bi-variate and multi-variate sensitivity analysis.

\item Function \code{collin} uses as input the sensitivity functions and
estimates the "collinearity" index for all possible parameter sets.
This multivariate sensitivity estimate measures approximate linear dependence
(and parameter identifiability).

\item Function \code{sensRange} produces 'envelopes' around the sensitivity
variables, as a function of the sensitivity parameters. This function can also
 be used to run simple "what-if" scenarios

\item Function \code{modMCMC} runs a markov chain monte carlo (Bayesian analysis).
It implements the delayed rejection - adaptive metropolis (DRAM) algorithm.

\item FME also contains unctions to generate multiple parameter values arranged
according to a grid (\code{Grid}), latin hypercube (\code{Latinhyper}),
multinormal (\code{Normal}) or uniform (\code{Uniform}) design.
}
}
\details{
  \tabular{ll}{
    Package: \tab FME\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2008-11-01\cr
    License: \tab GNU Public License 2 or above\cr
  }

}
\author{
  Karline Soetaert
}
%\references{
%% to be filled
%}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(modCost)
example(sensFun)


## open the directory with documents
browseURL(paste(system.file(package="FME"), "/doc", sep=""))

## the vignette
vignette("FME")

edit(vignette("FME"))
}
}

\keyword{ package }