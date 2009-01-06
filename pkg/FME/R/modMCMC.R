modMCMC <- function (f,p,...,
                     jump=NULL,
                     lower=-Inf, upper= +Inf,
                     prior=NULL,    # -2log probability or ssr  of parameters
                     var0=NULL,  # prior error variance
                     n0=NULL,      # prior accuracy for s2
                     niter=1000,
                     outputlength = niter,
                     updatecov = niter,
                     covscale = 2.4^2/length(p),
                     burninlength=0,
                     verbose=TRUE)

{
#======================
# 1. Initialisation...
#======================

# check input
 np  <- length(p)
 pnames <- names(p)

 if (updatecov == 0) updatecov <- niter
 isR       <- updatecov < niter  # will update the jump covariances
 nupdate   <- updatecov

# NewPars, function to generate new parameters or standarddeviation of pars
#---------

 NewParsMN <- function(p)   # multidimensional normal distribution...
   { pp<- p + rnorm(np)%*%R
     names(pp) <- pnames
     pp
    }
 if (is.null(jump))
    NewPars <- function(p) p + rnorm(np,sd=0.1)   else
 if (is.matrix (jump)) {
    R       <- chol(jump)
    NewPars <- NewParsMN
    nupdate <- updatecov
    }                                             else
 if (is.numeric(jump))
    NewPars <- function(p) p + rnorm(np,sd=jump)  else
                       NewPars <- jump
# Prior, function that returns the prior parameter probability
#---------
  if (is.null(prior)) Prior <- function(p) return(0) else
                      Prior <- function(p) return(prior(p))

# Func, function that returns the function value, either -2*log Probability or sum squares.
#---------
  uselogP <- is.null(var0)
  Pars    <- p
  limits  <- (lower != -Inf || upper != Inf)


  useSigma <- !is.null(var0)
  if (useSigma) divsigma <- 1/var0 else divsigma <- 1
  updateSigma <- (useSigma & !is.null(n0))
  
  FF      <- f(p,...) # First function call already done...
  if (is.numeric(FF)) N <- length(FF) else N <- nrow(FF$residual)
  PP      <- Prior(p)
  
  
  if (is.numeric(FF) & !useSigma)
     Func <- function(p,...)
     {
     PP   <<- Prior(p)
     FF   <<-f(p,...)
     return(0.5*(FF + PP))  #-log(probability)
     }     else

  if (!useSigma)      # use -2*log(probability)
     Func <- function(p,...)
     {
     PP   <<- Prior(p)
     FF   <<- f(p,...)$minlogp
     return(0.5*(FF + PP))  #-log(probability)
     }    else

     Func <- function(p,...)
     {
     PP       <<- Prior(p)
     FF       <<- f(p,...)  # use sum of squared residuals as -2 log(prob)...
     if(is.numeric(FF))      FF <<-sum(FF^2) else FF <<- FF$var$SSR.unweighted
     if(updateSigma)   divsigma <<- rgamma(1,shape=0.5*(n0+N),rate=0.5*(n0*var0+FF))
     return(0.5*(sum(FF*divsigma) + PP))
     }

  if (is.numeric(FF)& !useSigma) Fini <- 0.5*(FF + PP)                 else
  if (!useSigma)                 Fini <- 0.5*(FF$minlogp + PP)         else
  if (is.numeric(FF))          Fini <- 0.5*(sum(FF^2)*divsigma + PP) else
                               Fini <- 0.5*(FF$model*divsigma + PP)
  PPold <- PP
  if (!useSigma) {if(is.numeric(FF))   FFold <- FF   else FFold <- FF$minlogp} else
  if (is.numeric(FF)) FFold <- sum(FF^2) else FFold <- FF$model

# Adapt function calls if limits are imposed.
#---------
  if (!limits)
  {                           # no need to change parameters
   Fun     <- function(p,...) Func(p,...)

   } else  {

   Lower <- rep(lower,len=length(p))
   Upper <- rep(upper,len=length(p))
   if(any (Lower >= p) || any(Upper <= p)) stop("cannot proceed: initial parameters not within ranges")
     # lower and upper bounds...
    lu   <- which(is.finite(Lower) & is.finite(Upper))
    Pars[lu] <- tan(pi*((p[lu]-Lower[lu])/(Upper[lu]-Lower[lu]) - 0.5))

   # just lower bounds...
    l <- which(is.finite(Lower) & !is.finite(Upper))
    Pars[l] <-log(p[l]-Lower[l])

   # just upper bounds...
    u <- which(!is.finite(Lower) & is.finite(Upper))
    Pars[u] <- log(-p[u]+Upper[u])

    Fun    <- function(p,...)                # param transformation in function call
    {
     PP     <- p
     PP[lu] <- Lower[lu]+(Upper[lu]-Lower[lu])*(atan(p[lu])/pi + 0.5)
     PP[l]  <- Lower[l]+exp(p[l])
     PP[u]  <- Upper[u]-exp(p[u])
     names(PP) <-pnames                #nlminb 'forgets' parameter names...
     Func(PP,...)
    }
   }

# The acceptance function ...
  Accept <- function(fnew,fold)
  {
  if (fnew < bestfunp)
       {bestfunp<<-funpnew ; bestPar<<-parnew}
  if (!updateSigma)
  {
    if (fnew < fold)
    {accept <- TRUE
     }else
     accept <- (runif(1) < exp(fold-fnew))
  } else {    # ! updateSigma
     tst    <- exp(-0.5*sum((FF-FFold)*divsigma) + 0.5*(PP-PPold))
     if (tst<0) accept<-FALSE else if (tst>1) accept <- TRUE else accept <- (runif(1) < tst)
     if (accept) {FFold <- FF; PPold <- PP}
     }
     return(accept)
  }
  
# MCMC jumps...
#======================
 ou <- ceiling((niter - burninlength)/outputlength)
 niter <- niter - (niter - burninlength)%%ou
 outputlength <- (niter - burninlength)%/%ou
 ou1 <- burninlength + ou

 pars      <- matrix(nc=outputlength,nr=np)     #arrays with results
 funpars   <- vector(length=outputlength)
 if (useSigma) sig       <- vector(length=outputlength)  else sig<-NULL

 naccepted <- 0
 pp        <- Pars
 funpp     <- Fini            # initial function value

 bestPar   <- pp
 bestfunp  <- funpp
 
 ii        <- 0
 iiacc     <- 0
 for ( i in 1:niter)
 {
  parnew  <- NewPars(pp)       # new parameter set
  funpnew <- Fun(parnew,...)   # probability of new parameter set

  if (is.infinite(funpnew)) accept <- FALSE else

  if(accept <- Accept(funpnew , funpp))
    {
    naccepted <- naccepted +1
    pp    <- parnew
    funpp <- funpnew
    }
  if (i == ou1)                # current parameter set saved...
    {
    ii<-ii+1
    if (accept) iiacc<- iiacc+1
    pars[,ii]  <-pp
    funpars[ii]<-funpp
    if (useSigma) sig[ii]<- 1/divsigma
    ou1 <- i + ou
    if (isR & ii > nupdate & iiacc > 100)   # update jump parameters
     {
       Covar <- (cov(t(pars[,1:ii]))  + diag(1e-16,np))  * covscale
#       CV <- rbind(CV,as.vector(Covar))
       R     <- chol(Covar)
       NewPars <- NewParsMN
       nupdate <- nupdate + updatecov
     }
    }
 }
 
 
#======================
  if(limits)  # backtransform
  {
    pars[lu,]<-Lower[lu]+(Upper[lu]-Lower[lu])*(atan(pars[lu,])/pi + 0.5)
    pars[l,] <-Lower[l]+exp(pars[l,])
    pars[u,] <-Upper[u]-exp(pars[u,])
    
    bestPar[lu]<-Lower[lu]+(Upper[lu]-Lower[lu])*(atan(bestPar[lu])/pi + 0.5)
    bestPar[l] <-Lower[l]+exp(bestPar[l])
    bestPar[u] <-Upper[u]-exp(bestPar[u])
  }
  if (is.null(pnames)) pnames <- paste("p",1:np,sep="")
  
  rownames(pars) <- pnames
  if (verbose)
       cat("number of accepted runs: ", naccepted,
            " out of ", niter, " (", 100 * naccepted/niter,
              "%) \n", sep = "")
  res <-list(pars=t(pars),funp=funpars,naccepted=naccepted,sig=sig,
             bestpar= bestPar,bestfunp=bestfunp)

  class(res) <- "modMCMC"
  return(res)
}

pairs.modMCMC <- function (x, gap = 0, upper.panel = NA, diag.panel = NA,
    ...)
{
    panel.cor <- function(x, y) text(x = mean(range(x)), y = mean(range(y)),
        labels = format(cor(x, y), digits = 2))
    panel.hist <- function(x, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 2))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$counts
        y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "cyan")
    }
    if (!is.null(upper.panel) && is.na(upper.panel))
        upper.panel <- panel.cor
    if (!is.null(diag.panel) && is.na(diag.panel))
        diag.panel <- panel.hist
    X <- x$pars
    labels <- colnames(X)
    pairs(X, diag.panel = diag.panel, labels = labels, gap = gap,
        upper.panel = upper.panel, ...)
}

plot.modMCMC <- function (x, Full=FALSE, ...)
{
  np <- NP <- ncol(x$pars)
  if (Full) np <- np +2
  nc <- ceiling(sqrt(np))
  nr <- ceiling(np/nc)
  
  mf <- par(mfrow=c(nr,nc))
  on.exit(par(mf))

  for(i in 1:NP) plot(x$pars[,i],type="l",main=colnames(x$pars)[i],
                      xlab="iter",ylab="",...)
  if (Full) {plot(x$funp,type="l",main="function",xlab="iter",ylab="")
             matplot(x$sig,type="l",main="model variances",xlab="iter",ylab="")
  }
}

summary.modMCMC <- function (object, ...)
{
  data.frame(rbind(
  mean=apply(object$pars,2,mean),
  sd  =apply(object$pars,2,sd),
  min =apply(object$pars,2,min),
  max =apply(object$pars,2,max),
  q025 =apply(object$pars,2,quantile,probs=0.25),
  q050 =apply(object$pars,2,quantile,probs=0.5),
  q075 =apply(object$pars,2,quantile,probs=0.75)
  ))
}

range.modMCMC <- function (MCMC,solver,num=100,sensvar=NULL,Full=TRUE,...)
{
  Srange <- sensRange(MCMC$pars,solver=solver,num=num,sensvar=sensvar,Full=Full,...)
  return(Srange)
}