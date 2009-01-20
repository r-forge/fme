modMCMC <- function (f,p,...,
                     jump=NULL,
                     lower=-Inf, upper= +Inf,
                     prior=NULL,    # -2log probability or ssr  of parameters
                     var0=NULL,     # initial error variance
                     n0=NULL,       # prior accuracy for initial error variance
                     niter=1000,
                     outputlength = niter,
                     burninlength=0,
                     updatecov = niter,
                     covscale = 2.4^2/length(p),
                     ntrydr=1,      # number of tries for delayed rejection procedure
                     drscale =NULL, # scale factor
                     verbose=TRUE)

{
#===============================================================================
# 1. Initialisation...
#===============================================================================
 setbound <- TRUE #toggles off parameter transformation... does not seem to work -
 # KARLINE:TEST and perhaps REMOVE
#-------------------------------------------------------------------------------
# check input
#-------------------------------------------------------------------------------
 npar   <- length(p)
 pnames <- names(p)

 if (length(lower)!= npar & length(lower)!=1)
   stop("length of 'lower' should be either 1 or equal to number of parameters")
 if (length(upper)!= npar & length(upper)!=1)
   stop("length of 'upper' should be either 1 or equal to number of parameters")

# iterations saved
 if (outputlength > (niter-burninlength)) outputlength <- niter-burninlength
 if (burninlength >= niter) stop("burninlenght is larger than niter")

# recalculate the number of iterations...
 ou <- ceiling((niter - burninlength)/outputlength)  # output interval
 niter <- niter - (niter - burninlength)%%ou
 outputlength <- (niter - burninlength)%/%ou
 ou1 <- burninlength + ou                            # first output iteration nr

# delayed rejection input
 ntrydr <- max(1,ntrydr)
 if(is.null(drscale)) drscale=c(1/5,1/4,1/3)
 ldr <- length(drscale)
 if (ldr<ntrydr) drscale <- c(drscale, rep(drscale[ldr],ntrydr-ldr))

# adaptive metropolis input
 if (updatecov == 0) updatecov <- niter
 isR       <- updatecov < niter  # will update the jump covariances
 nupdate   <- updatecov          # first time update of covariance matrix

#-------------------------------------------------------------------------------
# Generating new parameter values
#-------------------------------------------------------------------------------

# function to be used if new parameters are generated with multiD normal distribution...
 NewParsMN <- function(p,R)
   { pp <- as.vector(p + rnorm(npar)%*%R) # R is cholesky decomposition of covariance mat
     names(pp) <- pnames
     return(pp)
    }
 if (is.null(jump)) {
    jump  <- 0.1*p # default is a sd of 0.1
    NewPars <- function(p,...) p + rnorm(npar,sd=jump)
    R <- diag(jump)
    }
 else if (is.matrix (jump)) {
    R       <- chol(jump)
    NewPars <- NewParsMN
    }
 else if (is.numeric(jump)){
    NewPars <- function(p,...) p + rnorm(npar,sd=jump)
    ii <- which (jump<0)
    jump[ii] <- 0.1*p[ii]
    R <- diag(nr=npar,jump)
    }
 else {
    NewPars <- jump             # jump is a function
    if (ntrydr > 1) stop ("cannot combine jump function with delayed rejection")
  }
#-------------------------------------------------------------------------------
# Prior, function that returns -2 n* log(prior parameter probability)
#-------------------------------------------------------------------------------
  if (is.null(prior)) Prior <- function(p) return(0) else   # uninformative
                      Prior <- function(p) return(prior(p))

  PPnew  <- Prior(p)

#-------------------------------------------------------------------------------
# check function call: function returns either
# -2*log Probability, the model residuals, or a list of class modFit.
#-------------------------------------------------------------------------------

  # if var0 has a value; used to scale the SSR
  useSigma    <- !is.null(var0)
  if (useSigma) divsigma <- 1/var0 else divsigma <- 1
  updateSigma <- (useSigma & !is.null(n0))
  if (useSigma) lenvar0 <- length(var0)
 # First function call ...
  SSnew  <- f(p,...)
  if (is.numeric(SSnew)& !useSigma) {
     if (length(SSnew)>1)
     stop("if 'var0' is NULL then function 'f' should either return the -2*log of the model probability, OR an instance of class 'modCost'")
     }
  else if (is.numeric(SSnew)& useSigma) {
    N <- length(SSnew)
    if (length(SSnew)==1)
     stop("if 'var0' has a value, then function 'f' should return the model residuals OR an instance of class 'modCost'")
    if (length(var0) != 1 && length(var0) != length(SSnew))
     stop("length of 'var0' should either 1 or = length of model residuals")
     }
  else if (class(SSnew) != "modCost")
    stop("function 'f' should either return the -2*log of the model probability OR an instance of class 'modCost'")
  else if (!useSigma)
     N <- nrow(SSnew$residuals)    # total number of data points
  else N <- SSnew$var$N            # number of data points per variable

  useRes <- FALSE
  if (class(SSnew) == "modCost" & useSigma)
  {
    if (length(var0) == nrow(SSnew$residuals))
      useRes <- TRUE
    else if (length(var0) != 1 & length(var0) != nrow(SSnew$var))
      stop("function 'f' is not compatible with length of var0")
  }


# note:SSnew and PPnew are globals...
  if (is.numeric(SSnew) & !useSigma) # f returns -2*log(probability(model))
     Func <- function(p,...)  {
       PPnew   <<- Prior(p)
       SSnew   <<- f(p,...)
       return(0.5*(SSnew + PPnew))      #-log(p_model*p_params)
     }
  else if (!useSigma)            # f returns a list of class modFit
     Func <- function(p,...)  {
       PPnew   <<- Prior(p)
       SSnew   <<- f(p,...)$minlogp # select -2*log(probability)
       return(0.5*(SSnew + PPnew))     #-log(p_model*p_params)
     }
  else if (is.numeric(SSnew))       # sigma + use sum of squared residuals
     Func <- function(p,...)  {
       PPnew  <<- Prior(p)
       SSnew  <<- f(p,...)
       SSnew  <<- SSnew^2
       if(updateSigma)           # draw from gamma distribution; var0 is prior model variance
          divsigma <<- rgamma(lenvar0,shape=0.5*(n0+N),rate=0.5*(n0*var0+SSnew))
       return(0.5*(sum(SSnew*divsigma) + PPnew))
     }
  else  if (useRes)              # sigma and residuals
     Func <- function(p,...)  {
       PPnew  <<- Prior(p)
       SSnew  <<- f(p,...)
       SSnew  <<- (SSnew$residuals$res.unweighted)^2
       if(updateSigma)           # draw from gamma distribution; var0 is prior model variance
          divsigma <<- rgamma(lenvar0,shape=0.5*(n0+N),rate=0.5*(n0*var0+SSnew))
       return(0.5*(sum(SSnew*divsigma) + PPnew))
     }
  else                           # sigma and varresiduals
     Func <- function(p,...)  {
       PPnew  <<- Prior(p)
       SSnew  <<- f(p,...)
       SSnew  <<- SSnew$var$SSR.unweighted
       if(updateSigma)           # draw from gamma distribution; var0 is prior model variance
          divsigma <<- rgamma(lenvar0,shape=0.5*(n0+N),rate=0.5*(n0*var0+SSnew))
       return(0.5*(sum(SSnew*divsigma) + PPnew))
     }

# initial function values
  if (is.numeric(SSnew)& !useSigma)
       SSold <- SSnew
  else if (!useSigma)
       SSold <- SSnew$minlogp
  else if (is.numeric(SSnew))
       SSold <- SSnew^2
  else if (useRes)
       SSold <- (SSnew$residuals$res.unweighted)^2
  else SSold <- SSnew$var$SSR.unweighted

  PPold   <- PPnew
  sigold  <- divsigma
  funpold <- 0.5* (sum(SSold*divsigma) + PPold)
  
#-------------------------------------------------------------------------------
# Adapt function calls if limits are imposed.
#-------------------------------------------------------------------------------
  parold      <- p
  lower[is.na(lower)] <- -Inf
  upper[is.na(upper)] <- Inf
  limits   <- (any(lower != -Inf) || any(upper != Inf))

  if (setbound)
    Fun     <- function(p,...) Func(p,...)
  else if (!limits)                # no need to change parameters
    Fun     <- function(p,...) Func(p,...)
  else  {
    Lower <- rep(lower,len=npar)
    Upper <- rep(upper,len=npar)
    if(any (Lower >= p) || any(Upper <= p))
      stop("cannot proceed: initial parameters not within ranges")
    # lower and upper bounds...
    lu   <- which(is.finite(Lower) & is.finite(Upper))
    parold[lu] <- tan(pi*((p[lu]-Lower[lu])/(Upper[lu]-Lower[lu]) - 0.5))

    # just lower bounds...
    l <- which(is.finite(Lower) & !is.finite(Upper))
    parold[l] <-log(p[l]-Lower[l])

    # just upper bounds...
    u <- which(!is.finite(Lower) & is.finite(Upper))
    parold[u] <- log(-p[u]+Upper[u])

    Fun    <- function(p,...) {        # param transformation in function call
      P     <- p
      P[lu] <- Lower[lu]+(Upper[lu]-Lower[lu])*(atan(p[lu])/pi + 0.5)
      P[l]  <- Lower[l]+exp(p[l])
      P[u]  <- Upper[u]-exp(p[u])
      names(P) <-pnames                #nlminb 'forgets' parameter names...
      return(Func(P,...))
    }
   }

#-------------------------------------------------------------------------------
# The acceptance function ...
#-------------------------------------------------------------------------------
  Test <- function(fnew,fold)
  {
    if (fnew < bestfunp)  # new best value
       {bestfunp<<-funpnew ; bestPar<<-parnew}
    if (!useSigma)  {
       alpha <- exp(fold-fnew)
    }
    else {    # uses SSR, scaled by sigma
     alpha  <- exp(-0.5*sum((SSnew-SSold)*divsigma) - 0.5*(PPnew-PPold))
    }
    return(alpha)
  }

  Accept <- function(tst)
  {
     if (tst<0)
       accept<-FALSE
     else if (tst>1)
       accept <- TRUE
     else accept <- (runif(1) < tst)

     if (accept) {
       SSold <<- SSnew
       PPold <<- PPnew
       if (useSigma) sigold <<- divsigma
     }
    return(accept)
  }
  
# in case delayed rejection...

A_count <- 0

AlphaFun<- function(arg) {

# recursive acceptance function for delayed rejection
# x.p, y1.p, ... contain the parameter value
# x.ss, y1.ss, ... the sum of squares
# x.a, y1.a, ... past alpha probabilities

  stage <- length(arg)-1      # the stage we are in: all the y's
  A_count <- A_count+1

# recursively compute past alphas
  a1 <- 1
  a2 <- 1

  for (k in 1:(stage-1)) {
#    a1 <- a1*(1-AlphaFun(arg[1:(k+1)]))  # already have these alphas
    a1 <- a1*(1-arg[(k+1)]$a)  # already have these alphas

    a2 <- a2*(1-AlphaFun(arg[seq(from=stage+1,to=stage+1-k,by=-1)]))
    if  (a2==0)  return (0)
  }
  ss2  <- arg[stage+1]$ss
  ss1  <- arg[1]$ss
  pri2 <- arg[stage+1]$pri
  pri1 <- arg[1]$pri

  y  <- -0.5*(sum((ss2-ss1)*divsigma) + (pri2-pri1))

  for (k in 1:stage) {
  y <- y + qfun(k,arg)
 }
 y <- min(1, exp(y)*a2/a1)
 
return(y)
}

# gaussian n-th stage log proposal ratio
# log of q_i(yn,... yn-j)/q_i(x,y1,...yj)
qfun <- function(iq,arg) {
   LL <- function(x) sum(x*x)

   stage <- length(arg)-2
   if (stage==iq)             # we are symmetric...
     return(0)
   else {
    iR <- invR[[iq]]
    y1 <- arg[[1]]$p            #y1
    y2 <- arg[[iq+1]]$p         #yi
    y3 <- arg[[stage+1]]$p      #yn
    y4 <- arg[[stage-iq+1]]$p   #yn-i
    z <- -0.5*(LL((y4-y3)%*%iR)-LL((y2-y1)%*%iR))
   }
  return(z)
}
  
#===============================================================================
# 2. MCMC jumps...
#===============================================================================

# matrices/vectors with results
 pars      <- matrix(nc=outputlength,nr=npar)
 funpars   <- vector(length=outputlength)

 if (useSigma)
   sig     <- matrix(nr=outputlength,nc=lenvar0)
 else sig<-NULL

# best function and parameter values will be kept
 bestPar   <- parold
 bestfunp  <- funpold

# keep track
 naccepted <- 0               # number of saved parameters
 ii        <- 0               # counter to saved parameters
 naccsave  <- 0               # number of saved accepted parameters
 ipos      <- 0

# need to save initial runs during burnin, if adaptive metropolis algorithm
 if (isR & burninlength>0) SaveIni <- TRUE else SaveIni <- FALSE
 if (SaveIni) # during burnin, the proposal covariance matrix needs to be updated...
 {
  naccsave2 <- 0
  icov      <- 0
  ou1b      <- updatecov
 }

 funpnew <- vector(length=ntrydr)
# if delayed rejection
 if (ntrydr>1){
     Rw   <- list()
     invR <- list()
     trypath <- list()
     x <- list()
     y <- list()
     z <- list()
 }

 Rold <- 0
 Rnew <- 1

 for ( i in 1:niter)
 {
   accept <- FALSE

   parnew     <- NewPars(parold,R)      # new parameter set

   if (setbound & (any(parnew < lower) | any(parnew > upper)))
     accept <- FALSE
   else {
      funpnew[1] <- Fun(parnew,...)        # probability of new parameter set

      if (is.infinite(funpnew[1])){
        alpha  <- 0
        accept <- FALSE
        }
      else  {
         alpha  <- Test(funpnew[1] , funpold)
         accept <- Accept(alpha)
      }
   }
   itry   <- 1

#-------------------------------------------------------------------------------
#  delayed rejection
#-------------------------------------------------------------------------------

   if (ntrydr > 1 & !accept) {
     x$p   <- parold
     x$ss  <- SSold
     x$pri <- PPold
     x$a   <- 0

     y$p   <- parnew
     y$ss  <- SSnew
     y$pri <- PPnew
     y$a   <- alpha

     trypath[[1]] <- x
     trypath[[2]] <- y

    # Create new R's and the inverse of R, scaled - ONLY IF R HAS CHANGED
     if (Rnew != Rold){
       invR[[1]] <- solve(R)
       Rw[[1]]   <- R
       for (i in 2:ntrydr) {
           Rw[[i]]   <- Rw[[i-1]] * drscale[i-1]
           invR[[i]] <- solve(Rw[[i]])
       }
       Rold <- Rnew
     }

     while (! accept & itry < ntrydr) {
       itry <- itry + 1

       parnew  <- NewPars(parold,Rw[[itry]])
       if (setbound & (any(parnew < lower) | any(parnew > upper))) {
           accept <- FALSE
           z$a <- 0
           z$p <- parnew
           z$pri<- 0
           z$ss <- Inf
       } else {
         funpnew[itry] <- Fun(parnew,...)  # probability of new parameter set
         z$p   <- parnew
         z$ss  <- SSnew
         z$pri <- PPnew
       }
       trypath[[itry+1]] <- z

       if (is.infinite(funpnew[itry]))
       { accept <- FALSE
         break
       }
       else {
         alpha <- AlphaFun(trypath[1:(itry+1)])
         trypath[itry+1]$a <- alpha
         accept <- Accept(alpha)
       }
     }
  }
#-------------------------------------------------------------------------------
# saving if accepted
#-------------------------------------------------------------------------------
   if (accept) {
     naccepted <- naccepted +1
     parold    <- parnew           # new parameter set replaces the old one...
     funpold   <- funpnew[itry]
   }
# -----  burnin: save values only if covariance needs updating  ------
   if (SaveIni) {  # use burnin to update covariance...
       ipos <- ipos+1
       if (ipos > outputlength) {
          ipos <- 1
          icov <- outputlength
       }
       pars[,ipos]  <-parold

       if (accept) naccsave2 <- naccsave2+1

       if (i> ou1b & naccsave2>5 ){  #  5 is arbitrary...
         jj <- max(ipos,icov)
         Covar <- (cov(t(pars[,1:jj]))  + diag(1e-16,npar)) * covscale
         RR    <- try(chol(Covar))
         if (is.numeric(RR)){
           R <- RR
           NewPars <- NewParsMN   # this function becomes the default update function
           ou1b <- i + updatecov
           Rnew <- Rnew + 1
         }
       }
    }
# -----  past burnin  ------

  # current parameter set saved every nupdate steps ...
  if (i == ou1)  {

    ii         <-ii+1
    pars[,ii]  <-parold
    funpars[ii]<-funpold
    if (useSigma)
      sig[ii,]<- 1/sigold

    if (accept) naccsave<- naccsave+1
    ou1 <- i + ou
    # update jump parameters by estimating param covariance..
    if (isR & ii > nupdate & naccsave > 5) {
       SaveIni <- FALSE
       Covar <- (cov(t(pars[,1:ii]))  + diag(1e-16,npar)) * covscale
       RR     <- try(chol(Covar))

       if (is.numeric(RR)){
         R <- RR
         NewPars <- NewParsMN   # this function becomes the default update function
         nupdate <- ii + updatecov
       Rnew <- Rnew + 1       }
     }
   }
 }

 
#===============================================================================
# 3. finalisation
#===============================================================================
  if(limits & ! setbound)  # backtransform
  {
    pars[lu,]<-Lower[lu]+(Upper[lu]-Lower[lu])*(atan(pars[lu,])/pi + 0.5)
    pars[l,] <-Lower[l]+exp(pars[l,])
    pars[u,] <-Upper[u]-exp(pars[u,])
    
    bestPar[lu]<-Lower[lu]+(Upper[lu]-Lower[lu])*(atan(bestPar[lu])/pi + 0.5)
    bestPar[l] <-Lower[l]+exp(bestPar[l])
    bestPar[u] <-Upper[u]-exp(bestPar[u])
  }
  if (is.null(pnames)) pnames <- paste("p",1:npar,sep="")
  if (useSigma) {
    if (length(var0) == ncol(sig))
       colnames(sig) <- names(var0)
    else colnames(sig) <- paste("var",1:ncol(sig),del="")
  }

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



################################################################################

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
    if (Full) X <- cbind(X,fvalue=x$funp)
    labels <- colnames(X)
    pairs(X, diag.panel = panel.hist, labels = labels, gap = 0,
        lower.panel = panel.cor, ...)
}

plot.modMCMC <- function (x, Full=FALSE, ...)
{
  np <- NP <- ncol(x$pars)
  if (Full) np <- np +1
  if (Full & !is.null(x$sig)) np <- np +ncol(x$sig)
  nc <- ceiling(sqrt(np))
  nr <- ceiling(np/nc)
  
  mf <- par(mfrow=c(nr,nc))
  on.exit(par(mf))

  for(i in 1:NP) plot(x$pars[,i],type="l",main=colnames(x$pars)[i],
                      xlab="iter",ylab="",...)
  if (Full) plot(x$funp,type="l",main="function value",xlab="iter",ylab="")
  if (Full & !is.null(x$sig)) for ( i in 1:ncol(x$sig))
       plot(x$sig[,i],type="l",main=colnames(x$sig)[i],
                      xlab="iter",ylab="variance",log="y")

}

hist.modMCMC <- function (x, Full=FALSE, breaks=100,...)
{
  np <- NP <- ncol(x$pars)
  if (Full) np <- np +1
  if (Full & !is.null(x$sig)) np <- np +ncol(x$sig)
  nc <- ceiling(sqrt(np))
  nr <- ceiling(np/nc)

  mf <- par(mfrow=c(nr,nc))
  on.exit(par(mf))

  for(i in 1:NP) hist(x$pars[,i],main=colnames(x$pars)[i],
                      xlab="",ylab="-",breaks=breaks,freq=FALSE,...)
  if (Full) hist(x$funp,main="function value",xlab="",ylab="-",breaks=breaks,freq=FALSE,...)
  if (Full & !is.null(x$sig)) for ( i in 1:ncol(x$sig))
       hist(sqrt(x$sig[,i]),main="error std posterior",
                      xlab="",ylab=colnames(x$sig)[i],
                      breaks=breaks,freq=FALSE,...)

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

