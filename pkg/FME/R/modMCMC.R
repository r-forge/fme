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
#===============================================================================
# 1. Initialisation...
#===============================================================================

# check input
 np     <- length(p)
 pnames <- names(p)

 if (updatecov == 0) updatecov <- niter
 isR       <- updatecov < niter  # will update the jump covariances
 nupdate   <- updatecov          # first time update of covariance matrix

#-------------------------------------------------------------------------------
# NewPars, function to generate new parameters or standarddeviation of pars
#-------------------------------------------------------------------------------

 NewParsMN <- function(p)       # multidimensional normal distribution...
   { pp <- as.vector(p + rnorm(np)%*%R)    # R is cholesky decomposition of covariance mat
     names(pp) <- pnames
     return(pp)
    }
    
 if (is.null(jump))
    NewPars <- function(p) p + rnorm(np,sd=0.1)   else

 if (is.matrix (jump)) {
    R       <- chol(jump)
    NewPars <- NewParsMN
    }                                             else
 if (is.numeric(jump))
    NewPars <- function(p) p + rnorm(np,sd=jump)  else
    NewPars <- jump             # it is a function

#-------------------------------------------------------------------------------
# Prior, function that returns -2 n* log(prior parameter probability)
#-------------------------------------------------------------------------------
  if (is.null(prior)) Prior <- function(p) return(0) else   # uninformative
                      Prior <- function(p) return(prior(p))

#-------------------------------------------------------------------------------

  # if var0 has a value; used to scale the SSR
  useSigma <- !is.null(var0)
  if (useSigma) divsigma <- 1/var0 else divsigma <- 1
  updateSigma <- (useSigma & !is.null(n0))

 # First function call ...
  FF  <- f(p,...)
  if (is.numeric(FF))
     N <- length(FF)            # total number of data points
  else if (!useSigma)
     N <- nrow(FF$residuals)     # total number of data points
  else N <- FF$var$N            # number of data points per variable

  PP  <- Prior(p)
#-------------------------------------------------------------------------------
# Func, function that returns the function value,
# either -2*log Probability, the model residuals, or a list of class modFit.
# note:FF and PP are locals...
#-------------------------------------------------------------------------------
  if (is.numeric(FF) & !useSigma) # f returns -2*log(probability(model))
     Func <- function(p,...)  {
       PP   <<- Prior(p)
       FF   <<- f(p,...)
       return(0.5*(FF + PP))      #-log(p_model*p_params)
     }
  else if (!useSigma)            # f returns a list of class modFit
     Func <- function(p,...)  {
       PP   <<- Prior(p)
       FF   <<- f(p,...)$minlogp # select -2*log(probability)
       return(0.5*(FF + PP))     #-log(p_model*p_params)
     }
  else
     Func <- function(p,...)  {
       PP  <<- Prior(p)
       FF  <<- f(p,...)
       if(is.numeric(FF))        # use sum of squared residuals
          FF <<-sum(FF^2)
       else FF <<- FF$var$SSR.unweighted
       if(updateSigma)           # draw from gamma distribution; var0 is prior model variance
          divsigma <<- rgamma(1,shape=0.5*(n0+N),rate=0.5*(n0*var0+FF))
       return(0.5*(sum(FF*divsigma) + PP))
     }

# initial function values
  if (is.numeric(FF)& !useSigma)
       FFold <- FF
  else if (!useSigma)
       FFold <- FF$minlogp
  else if (is.numeric(FF))
       FFold <- sum(FF^2)
  else FFold <- FF$var$SSR.unweighted

  PPold <- PP
  funpp <- 0.5* (sum(FFold*divsigma) + PPold)
  
#-------------------------------------------------------------------------------
# Adapt function calls if limits are imposed.
#-------------------------------------------------------------------------------
  pp      <- p
  limits   <- (lower != -Inf || upper != Inf)

  if (!limits)                # no need to change parameters
    Fun     <- function(p,...) Func(p,...)
  else  {
    Lower <- rep(lower,len=length(p))
    Upper <- rep(upper,len=length(p))
    if(any (Lower >= p) || any(Upper <= p))
      stop("cannot proceed: initial parameters not within ranges")
    # lower and upper bounds...
    lu   <- which(is.finite(Lower) & is.finite(Upper))
    pp[lu] <- tan(pi*((p[lu]-Lower[lu])/(Upper[lu]-Lower[lu]) - 0.5))

    # just lower bounds...
    l <- which(is.finite(Lower) & !is.finite(Upper))
    pp[l] <-log(p[l]-Lower[l])

    # just upper bounds...
    u <- which(!is.finite(Lower) & is.finite(Upper))
    pp[u] <- log(-p[u]+Upper[u])

    Fun    <- function(p,...) {        # param transformation in function call
      PP     <- p
      PP[lu] <- Lower[lu]+(Upper[lu]-Lower[lu])*(atan(p[lu])/pi + 0.5)
      PP[l]  <- Lower[l]+exp(p[l])
      PP[u]  <- Upper[u]-exp(p[u])
      names(PP) <-pnames                #nlminb 'forgets' parameter names...
      return(Func(PP,...))
    }
   }

#-------------------------------------------------------------------------------
# The acceptance function ...
#-------------------------------------------------------------------------------
  Accept <- function(fnew,fold)
  {
    if (fnew < bestfunp)  # new best value
       {bestfunp<<-funpnew ; bestPar<<-parnew}
    if (!updateSigma)  {
      if (fnew < fold)
        accept <- TRUE
      else
       accept <- (runif(1) < exp(fold-fnew))
    }
    else {    # updateSigma
     tst    <- exp(-0.5*sum((FF-FFold)*divsigma) + 0.5*(PP-PPold))
     if (tst<0)
       accept<-FALSE
     else if (tst>1)
       accept <- TRUE
     else accept <- (runif(1) < tst)

     if (accept) {
       FFold <- FF
       PPold <- PP
     }
    }
    return(accept)
  }
  
#===============================================================================
# MCMC jumps...
#===============================================================================
# recalculate the number of iterations...
 ou <- ceiling((niter - burninlength)/outputlength)  # output interval
 niter <- niter - (niter - burninlength)%%ou
 outputlength <- (niter - burninlength)%/%ou
 ou1 <- burninlength + ou                            # first output iteration nr

# matrices/vectors with results
 pars      <- matrix(nc=outputlength,nr=np)
 funpars   <- vector(length=outputlength)
 if (useSigma)
   sig     <- vector(length=outputlength)
 else sig<-NULL

 if (isR & burninlength>0) SaveIni <- TRUE else SaveIni <- FALSE

 bestPar   <- pp              # best function and parameter values will be kept
 bestfunp  <- funpp
 naccepted <- 0               # number of saved parameters
 ii        <- 0               # counter to saved parameters
 naccsave  <- 0               # number of saved accepted parameters
 ipos      <- 0

 if (SaveIni)
 {
  naccsave2 <- 0
  icov      <- 0
  ou1b      <- updatecov
 }
 for ( i in 1:niter)
 {
  parnew  <- NewPars(pp)      # new parameter set
  funpnew <- Fun(parnew,...)  # probability of new parameter set
  if (is.infinite(funpnew))
    accept <- FALSE
  else if(accept <- Accept(funpnew , funpp)) {
    naccepted <- naccepted +1
    pp    <- parnew           # new parameter set replaces the old one...
    funpp <- funpnew
  }

  if (SaveIni) {  # use burnin to update covariance...
       ipos <- ipos+1
       if (ipos > outputlength) {ipos <- 1; icov <- outputlength}
       pars[,ipos]  <-pp
       if (accept) naccsave2 <- naccsave2+1
       if (i> ou1b & naccsave2>5 ){
       jj <- max(ipos,icov)
       Covar <- (cov(t(pars[,1:jj]))  + diag(1e-16,np)) * covscale
       RR    <- try(chol(Covar))
       if (is.numeric(RR))
       { R <- RR
         NewPars <- NewParsMN   # this function becomes the default update function
         ou1b <- i + updatecov
       }
       }
  }

  # current parameter set saved...
  if (i == ou1)  {
    ii<-ii+1
    if (accept) naccsave<- naccsave+1
    pars[,ii]  <-pp
    funpars[ii]<-funpp
    if (useSigma) sig[ii]<- 1/divsigma
    ou1 <- i + ou
    # update jump parameters by estimating param covariance..
    if (isR & ii > nupdate & naccsave > 5) {
       SaveIni <- FALSE
       Covar <- (cov(t(pars[,1:ii]))  + diag(1e-16,np)) * covscale
       RR     <- try(chol(Covar))
       if (is.numeric(RR)){
         R <- RR
         NewPars <- NewParsMN   # this function becomes the default update function
         nupdate <- ii + updatecov
       }
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

pairs.modMCMC <- function (x, Full=FALSE, ...)
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
    X <- x$pars
    if (Full) X <- cbind(X,x$funp)
    labels <- colnames(X)
    pairs(X, diag.panel = panel.hist, labels = labels, gap = 0,
        upper.panel = panel.cor, ...)
}

plot.modMCMC <- function (x, Full=FALSE, ...)
{
  np <- NP <- ncol(x$pars)
  if (Full) np <- np +1
  if (Full & !is.null(x$sig)) np <- np +1
  nc <- ceiling(sqrt(np))
  nr <- ceiling(np/nc)
  
  mf <- par(mfrow=c(nr,nc))
  on.exit(par(mf))

  for(i in 1:NP) plot(x$pars[,i],type="l",main=colnames(x$pars)[i],
                      xlab="iter",ylab="",...)
  if (Full) plot(x$funp,type="l",main="function value",xlab="iter",ylab="")
  if (Full & !is.null(x$sig))  matplot(x$sig,type="l",main="model variances",xlab="iter",ylab="")

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

