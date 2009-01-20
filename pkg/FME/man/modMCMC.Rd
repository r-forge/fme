\name{modMCMC}
\alias{modMCMC}
\alias{summary.modMCMC}
\alias{plot.modMCMC}
\alias{pairs.modMCMC}
\alias{hist.modMCMC}
\title{Constrained markov chain monte carlo}
\description{performs a Markov Chain Monte Carlo simulation,
according to the adaptive metropolis (AM) algorithm and including a delayed rejection (DR) procedure}
\usage{modMCMC(f, p, ..., jump=NULL,  lower=-Inf, upper=+Inf,
prior=NULL, var0 = NULL, n0 = NULL,
niter=1000, outputlength = niter, burninlength=0,
updatecov=niter, covscale = 2.4^2/length(p),
ntrydr=1, drscale=NULL, verbose=TRUE)

\method{summary}{modMCMC}(object, \dots)
\method{pairs}{modMCMC}(x, Full=FALSE, \dots)
\method{hist}{modMCMC}(x, Full=FALSE, breaks=100, \dots)
\method{plot}{modMCMC}(x, Full=FALSE, \dots)
                     }
\arguments{
  \item{f }{A function to be minimized, with first argument the vector of parameters
            which should be varied,
            It should return either -2*log of the model probability (one value),
            the model residuals, or an element of class \emph{modCost} (as returned by a call to \code{\link{modCost}}.}
  \item{p }{Initial values for the parameters to be optimized over}
  \item{...}{additional arguments passed to function \code{f} (modMCMC) or to the methods}
  \item{jump }{jump length (for bounded parameters, in transformed space), either a \emph{number} or a
  \emph{vector} with length = total number of parameters, or a \emph{function} that takes as input the current
  values of the parameters and produces as output the perturbed parameters. See details}
  \item{prior }{prior probability of the parameters, either a function that is called as \code{prior(p)}
    or \code{NULL}; in the latter case a non-informative prior is used (i.e. all parameters are equally likely)}
  \item{var0 }{Initial model variance; if \code{NULL},
  it is assumed that initial model variance is 1, and the return element from \code{f} is -2*log (probability),
  else it is assumed that the return element from \code{f} is the sum of squared residuals. - see details}
  \item{n0 }{"accuracy" for the initial model variance - see details}
  \item{lower }{Lower bounds on the parameters; if unbounded set equal to \code{-Inf}}
  \item{upper }{Upper bounds on the parameters; if unbounded set equal to \code{Inf}}
  \item{niter }{number of iterations for MCMC}
  \item{outputlength }{number of iterations kept in the output; should be smaller than \code{niter}}
  \item{updatecov }{number of iterations after which the parameter covariance matrix is (re)evaluated based on
  the parameters kept thus far, and used to update the MCMCjumps}
  \item{covscale }{scale factor for the parameter covariance matrix, used to perform the MCMCjumps}
  \item{burninlength }{number of initial iterations to be removed from output}
  \item{ntrydr }{number of tries for the delayed rejection procedure}
  \item{drscale }{for each try during delayed rejection, the cholesky decomposition of the
  proposal matrix is scaled with this amount; if NULL, it is assumed to be c(0.2,0.25, 0.333, 0.333,...)}
  \item{verbose }{if TRUE: prints extra output}
  \item{object }{an object of class \code{modMCMC}}
  \item{x }{an object of class \code{modMCMC}}
  \item{Full }{If TRUE then not only the parameters will be plotted, but also the function value and (if appropriate) the model variance(s)}
  \item{breaks }{the number of cells for the histogram}
                      }
\value{
  a list of class \emph{modMCMC} containing the results as returned from the markov chain.
  
  This includes the following:
  
  \item{par }{an array with dimension (\code{outputlenght},length(p)), containing the parameters of the MCMC}
  \item{logp }{vector with the logarithm of the probability for each row in \code{par}}
  \item{naccepted }{the number of accepted runs}
  \item{bestpar }{the parameter set that gave the highest probability}
  \item{bestlogp }{the logarithm of the highest probability (corresponding to \code{bestpar}}

  The list returned by \code{modMCMC} has methods for the generic functions
   \code{\link{summary}},  \code{\link{plot}},   \code{\link{pairs}}- see note.

}
\details{
Note that arguments after ... must be matched exactly.

R-function \code{f} is called as \code{f(p,...)}. It should return either -2 times the log
likelihood of the model, one value, or an item of class
\code{modFit} (as created by function \code{\link{modFit}}.
In the latter case, the model SSR (which is equal to -2*log likelhood for Gaussian distributions)
or the variable SSR, will be used.

\code{var0} gives the variance of the measured components. It is either one value, or one value for each
observed variable, or one for each observed value.
When \code{var0} is not NULL, then \code{f} is assumed to return the model residuals OR an instance of
class \code{modCost}. When not NULL, then \code{f} should return either -2*log(probability of the model) an instance of
class \code{modCost}.

If \code{n0} is 0 (the default) then the variances are
assumed to be fixed. Typically, these variances can be estimated from the sum of squares of
fitted residuals. (e.g. as reported in \code{modFit}).

If \code{n0} is >0, then the variances are sampled as conjugate priors from the inverse
gamma distribution with parameters \code{var0} and \code{n0}. Larger values of \code{n0} keep
these samples closer to \code{var0}.

\code{modMCMC} implements the Metropolis-Hastings method. The proposal distribution, which
is used to generate new parameter values is the (multidimensional) Gaussian density distribution, with standard
deviation given by \code{jump}.

\code{jump} can be either one value, a vector of length = number of parameters or a
parameter covariance matrix (nrow=ncol=number parameters).

The jump parameter, \code{jump} thus determines how much the new parameter set will deviate from the old one.
It can be either one value, a vector or a matrix.

If \code{jump} is one value, or a vector, then the new parameter values are generated
by sampling a normal distribution with standard deviation equal to \code{jump}.
A larger value will lead to larger jumps in the parameter space, but acceptance of new points
can get very low. Smaller jump lengths increase the acceptance rate,
but the algorithm may move too slowly, and too many runs may be needed to scan the parameter space.

If \code{jump} is NULL, then the jump length is taken as 10\% of the parameter value.

\code{jump} can also be a proposal covariance matrix. In this case, the new parameter
values are generated by sampling a multidimensional normal distribution.

Finally, \code{jump} can also be an R-function that takes as input the current values
of the parameters and returns the new parameter values.

Two methods are implemented to increase the number of accepted runs.

1. In the "adaptive metropolis" method, new parameters are generated with a
covariance matrix that is estimated from the parameters generated (and saved) thus far.
Setting \code{updatecov} smaller to \code{niter} will trigger this function.
In this case, every \code{updatecov} iterations, the jump covariance matrix will be estimated
from the covariance matrix of the saved parameter values. The covariance matrix is scaled
with $(2.4^2/npar)$ where npar is the number of parameters, unless \code{covscal} has
been given a different value.
Thus, \eqn{Jump= ( cov(p1,p2,....pn) \cdot diag(np,+1e^{-16}) /(2.4^2/npar)} where the small number $1e^{-16}$ is added
on the diagonal of the covariance matrix to prevent it from becoming singular.

When using the adaptive Metropolis method, it is best to start with a small value of the
jump length.

2. In the "delayed rejection" method, new parameter values are tried upon rejection.
The number of tries is set with \code{ntry}. Setting this equal to 1 toggles off delayed rejection.
New parameters are generated by scaling the jump covariance matrix with a factor as specified in
\code{drscale}.

Convergence can be checked via \code{plot}, which plots the values of all parameters.
There should be no visible drift.

It is possible to select one jump length, or one for each parameter, or to define a function
that takes as input the current parameters and that produces the new parameter set.

In the first two cases, new parameter values will be generated by drawing from a normal distribution
with mean 0 and standard deviation = \code{jump}.

The \code{burninlength} is the number of initial steps that is not included in the output.
It can be useful if the initial value of the parameters is far from the optimal value.
Starting the MCMC with the best fit parameter set will alleviate the need for \code{burninlength}.




}

\note{
The following methods are provided:
\itemize{
\item{summary }{produces summary statistics of the MCMC results}
\item{plot }{plots the MCMC results, for all parameters. Use it to check convergence}
\item{pairs }{produces a pairs plot of the MCMC results; overrides the default
  \code{gap = 0}, \code{upper.panel = NA}, and \code{diag.panel}}
}
}

\seealso{\code{\link{modFit}} for constrained model fitting}
\examples{

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# testing: sampling a distribution
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# test 1: sampling a normal distribution, mean=1:3, sd = 0.25

NN <- function(p)
{
mu <- c(1,2,3)
-2*sum(log(dnorm(p,mean=mu,sd=0.1)))      #-2*log(probability)
}

# no bounds imposed...

MCMC <- modMCMC (f=NN,p=0:2, niter=10000,
                 outputlength=1000,jump=0.5)

# More accepted values by updating the jump covariance matrix...
MCMC <- modMCMC (f=NN,p=0:2, niter=10000,updatecov=100,
                 outputlength=1000,jump=0.1)
summary(MCMC)

plot(MCMC)   # check convergence
pairs(MCMC)

# bounds imposed...
NN <- function(p)
{
mu <- c(1,2,3)
-2*sum(log(dnorm(p,mean=mu,sd=0.5)))      #-2*log(probability)
}
MCMC2 <- modMCMC (f=NN,p=0:2+0.1, niter=10000,burninlength=2500,
                  outputlength=1000,jump=0.7, lower=-3:-1,upper=6:8)
plot(MCMC2)
summary(MCMC2)

# test 2: sampling a log-normal distribution, log mean=1:4, log sd = 1

NL <- function(p)
{
mu <- 1:4
-2*sum(log(dlnorm(p,mean=mu,sd=1)))      #-2*log(probability)
}
MCMCl <- modMCMC (f=NL,p=log(1:4),niter=10000,outputlength=1000,jump=5)
plot(MCMCl)   # bad convergence

MCMCl <- modMCMC (f=NL,p=log(1:4),niter=10000,outputlength=1000,jump=2^(2:5))
plot(MCMCl)   # better convergence but CHECK it!
pairs(MCMCl)
colMeans(log(MCMCl$pars))
sd(log(MCMCl$pars))

MCMCl <- modMCMC (f=NL,p=rep(1,4),niter=5000,outputlength=1000,
         jump=5,updatecov=100)
plot(MCMCl)
colMeans(log(MCMCl$pars))
sd(log(MCMCl$pars))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# example : Fitting a Monod (Michaelis-Menten) function to data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# the observations
#---------------------
Obs <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour
plot(Obs,pch=16,cex=2,xlim=c(0,400),ylim=c(0,0.15),
     xlab="mg COD/l",ylab="1/hr",main="Monod")

# the Monod model
#---------------------
Model <- function(p,x)   data.frame(x=x,N=p[1]*x/(x+p[2]))

# Fitting the model to the data
#---------------------
# define the residual function
Residuals  <- function(p) (Obs$y-Model(p,Obs$x)$N)  #... model residuals

# use modFit to find parameters
P      <- modFit(f=Residuals,p=c(0.1,1))

# plot best-fit model
x      <-0:375
lines(Model(P$par,x))

# summary of fit
sP    <- summary(P)
sP[]
print(sP)

# Running an MCMC
#---------------------
# estimate parameter covariances
# (to efficiently generate new parameter values)
Covar   <- sP$cov.scaled * 2.4^2/2

# the model variance
s2prior <- sP$modVariance

# set nprior = 0 to avoid updating model variance
MCMC <- modMCMC(f=Residuals,p=P$par,jump=Covar,niter=10000,
                var0=s2prior,n0=NULL,updatecov=100)

plot(MCMC,Full=TRUE)
pairs(MCMC)

cor(MCMC$pars)

cov(MCMC$pars)   # covariances by MCMC
sP$cov.scaled    # covariances by hessian of fit

x  <- 1:400
SR <- summary(sensRange(parInput=MCMC$pars,func=Model,x=x))
plot(SR, xlab="mg COD/l",ylab="1/hr",main="Monod")
points(Obs,pch=16,cex=1.5)

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\references{
H. Haario, E. Saksman and J. Tamminen, 2001. An adaptive Metropolis algorithm Bernoulli 7, pp. 223-242.

H. Haario, E. Saksman and J. Tamminen, 2005. Componentwise adaptation for high dimensional MCMC.
Computational Statistics 20(2), 265-274.

A. Gelman, JB Varlin, HS Stern and DB Rubin, 2004. Bayesian Data Analysis, second edition.
  Chapman and Hall,Boca Raton.
}
\keyword{utilities}

