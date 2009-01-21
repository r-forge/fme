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
                     ntrydr=1,      # maximal number of tries for delayed rejection procedure
                     drscale =NULL, # scale factor in DR
                     verbose=TRUE)
# Part of this R-code (the DR-part) is based on the matlab code (Marco Laine)
# as available from http://www.helsinki.fi/~mjlaine/mcmc/
{
#===============================================================================
# 1. Initialisation...
#===============================================================================

#-------------------------------------------------------------------------------
# check input
#-------------------------------------------------------------------------------
 npar   <- length(p)
 pnames <- names(p)
 
 if (length(lower)!= npar & length(lower)!=1)
   stop("length of 'lower' should be either 1 or equal to number of parameters")
 if (length(upper)!= npar & length(upper)!=1)
   stop("length of 'upper' should be either 1 or equal to number of parameters")

  lower[is.na(lower)] <- -Inf
  upper[is.na(upper)] <- Inf
  limits   <- (any(lower != -Inf) || any(upper != Inf))

# number of iterations saved
 if (outputlength > (niter-burninlength)) outputlength <- niter-burninlength

# recalculate the number of iterations ~ burnin and output...
 if (burninlength >= niter) stop("burninlenght is larger than niter")
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
 isR       <- updatecov < niter  # if TRUE: will update the jump covariances
 nupdate   <- updatecov          # first time update of covariance matrix

#-------------------------------------------------------------------------------
# Generating new parameter values
#-------------------------------------------------------------------------------
# function to be used if new parameters are generated with multiD normal distribution...
# (this may not be the case at first, but come into play if updatecov < niter
 NewParsMN <- function(p,R)
   { pp <- as.vector(p + rnorm(npar)%*%R) # R is cholesky decomposition of cov mat
     names(pp) <- pnames
     return(pp)
    }

# Check jump input (to generate new parameters)
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
# model variance: initial value and accurracy...
#-------------------------------------------------------------------------------
# if var0 has a value, it is used to scale the SSR; divsigma = 1/var
  useSigma    <- !is.null(var0)
  if (useSigma) divsigma <- 1/var0 else divsigma <- 1

# if n0 has a value, then sigma will be updated (updateSigma will be TRUE)
  updateSigma <- (useSigma & !is.null(n0))

  lenvar0 <- length(var0)

#-------------------------------------------------------------------------------
# check function call: function returns either
# -2*log Probability, the model residuals, or a list of class modFit.
#-------------------------------------------------------------------------------
 # First function call ...
  SSnew  <- f(p,...)
  if (is.numeric(SSnew)& !useSigma) {
     if (length(SSnew)>1)
     stop("if 'var0' is NULL then function 'f' should either return -2*log of the model probability, OR an instance of class 'modCost'")
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

#-------------------------------------------------------------------------------
# The function call used in the MCMC
#-------------------------------------------------------------------------------

# note:SSnew and PPnew are globals...
  if (is.numeric(SSnew) & !useSigma)   # f returns -2*log(probability(model))
     Fun <- function(p,...)  {
       PPnew   <<- Prior(p)
       SSnew   <<- f(p,...)
       return(0.5*(SSnew + PPnew))     #-log(p_model*p_params)
     }
  else if (!useSigma)                  # f returns a list of class modFit
     Fun <- function(p,...)  {
       PPnew   <<- Prior(p)
       SSnew   <<- f(p,...)$minlogp    # select -2*log(probability)
       return(0.5*(SSnew + PPnew))     #-log(p_model*p_params)
     }
  else if (is.numeric(SSnew))          # sigma + use sum of squared residuals
     Fun <- function(p,...)  {
       PPnew  <<- Prior(p)
       SSnew  <<- f(p,...)
       SSnew  <<- SSnew^2
       if(updateSigma)                 # draw from gamma distribution; var0 is prior model variance
          divsigma <<- rgamma(lenvar0,shape=0.5*(n0+N),rate=0.5*(n0*var0+SSnew))
       return(0.5*(sum(SSnew*divsigma) + PPnew))
     }
  else  if (useRes)              # sigma and residuals
     Fun <- function(p,...)  {
       PPnew  <<- Prior(p)
       SSnew  <<- f(p,...)
       SSnew  <<- (SSnew$residuals$res.unweighted)^2
       if(updateSigma)           # draw from gamma distribution; var0 is prior model variance
          divsigma <<- rgamma(lenvar0,shape=0.5*(n0+N),rate=0.5*(n0*var0+SSnew))
       return(0.5*(sum(SSnew*divsigma) + PPnew))
     }
  else                           # sigma and varresiduals
     Fun <- function(p,...)  {
       PPnew  <<- Prior(p)
       SSnew  <<- f(p,...)
       SSnew  <<- SSnew$var$SSR.unweighted
       if(updateSigma)           # draw from gamma distribution; var0 is prior model variance
          divsigma <<- rgamma(lenvar0,shape=0.5*(n0+N),rate=0.5*(n0*var0+SSnew))
       return(0.5*(sum(SSnew*divsigma) + PPnew))
     }

#-------------------------------------------------------------------------------
# initial values for f, SS, PP, sigma...
#-------------------------------------------------------------------------------
  if (is.numeric(SSnew)& !useSigma)
       SSold <- SSnew
  else if (!useSigma)
       SSold <- SSnew$minlogp
  else if (is.numeric(SSnew))
       SSold <- SSnew^2
  else if (useRes)
       SSold <- (SSnew$residuals$res.unweighted)^2
  else SSold <- SSnew$var$SSR.unweighted

  parold  <- p
  PPold   <- PPnew     # prior
  sigold  <- divsigma
  funpold <- 0.5* (sum(SSold*divsigma) + PPold)
  SSnew   <- Inf

#-------------------------------------------------------------------------------
# The acceptance function ...
#-------------------------------------------------------------------------------
  Test <- function(fnew,fold)
  {
    if (is.na(fnew) | is.infinite(fnew)) return(0)
    if (fnew < bestfunp)  # new best value - save it!
       {bestfunp<<-fnew ; bestPar<<-parnew}
    if (!useSigma)  {
       test <- exp(fold-fnew)
    }
    else {    # uses SSR, scaled by sigma
     test  <- exp(-0.5*sum((SSnew-SSold)*divsigma) - 0.5*(PPnew-PPold))
    }
    return(test)
  }

  Accept <- function(tst)
  {
     if (is.nan(tst)|is.na(tst))
       accept <- FALSE
     else if (tst<0)
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
  
#-------------------------------------------------------------------------------
# in case delayed rejection procedure ...
#-------------------------------------------------------------------------------

A_count  <- 0       # counter to number of alpha entrances
dr_steps <- 0       # counter to number of dr steps...


AlphaFun <- function(arg) {
#-------------------------------------------------------------------------------
# recursive acceptance function for delayed rejection
# x.p, y1.p, ... contain the parameter value
# x.ss, y1.ss, ... the sum of squares
#-------------------------------------------------------------------------------

  stage <- length(arg)-1      # the stage we are in
  A_count <<- A_count+1

# recursively compute past alphas
  a1 <- 1
  a2 <- 1

  if (stage > 1){
    for (k in 1:(stage-1)) {
    a1 <- a1*(1-AlphaFun(arg[1:(k+1)]))  # already have these alphas
#      a1 <- a1*(1-arg[(k+1)]$a)  # already have these alphas   k+1???

    a2 <- a2*(1-AlphaFun(arg[seq(from=stage+1,to=stage+1-k,by=-1)]))
#      a2 <- a2*(1-Aarg[(k+1)]$b) # these alphas have been calculated
      if  (is.null(a2) | is.na(a2))  return (0)
      if  (a2==0)  return (0)
    }
  }
  ss2  <- arg[[stage+1]]$ss
  ss1  <- arg[[1]]$ss
  pri2 <- arg[[stage+1]]$pri
  pri1 <- arg[[1]]$pri

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

   stage <- length(arg)-1   # CHANGED FROM MATLAB CODE: was -2
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
 pars      <- matrix(nc=outputlength,nr=npar)  # parameter values
 funpars   <- vector(length=outputlength)      # function value

 if (useSigma)
   sig     <- matrix(nr=outputlength,nc=lenvar0)  # model variances
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
 if (isR & burninlength>0)    # during burnin, the proposal covariance matrix needs to be updated...
 { SaveIni    <- TRUE
  naccsave2   <- 0
  icov        <- 0
  burnupdate  <- updatecov
 } else SaveIni <- FALSE

# if delayed rejection
 if (ntrydr>1){
     Rw   <- list()      # The cholesky decompositions used
     invR <- list()      # The inverses of the choleski decompositions
     trypath <- list()   # The parameter tries
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

   if (any(parnew < lower) | any(parnew > upper))
     accept <- FALSE
   else {
      funpnew <- Fun(parnew,...)        # probability of new parameter set

      if (is.infinite(funpnew)){
        alpha  <- 0
        accept <- FALSE
        }
      else  {
         alpha  <- Test(funpnew , funpold)
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

     y$p   <- parnew
     y$ss  <- SSnew
     y$pri <- PPnew

     trypath[[1]] <- x
     trypath[[2]] <- y

    # Create new R's and the inverse of R, scaled - ONLY IF R HAS CHANGED
     if (Rnew != Rold){
       invR[[1]] <- solve(R)
       Rw[[1]]   <- R
       for (j in 2:ntrydr) {
           Rw[[j]]   <- Rw[[j-1]] * drscale[j-1]
           invR[[j]] <- solve(Rw[[j]])
       }
       Rold <- Rnew
     }

     while (! accept & itry < ntrydr) {
       itry <- itry + 1

       parnew  <- NewPars(parold,Rw[[itry]])
       if (any(parnew < lower) | any(parnew > upper)) {
           accept <- FALSE
           z$p <- parnew
           z$pri<- 0
           z$ss <- Inf
       } else {
         funpnew <- Fun(parnew,...)  # probability of new parameter set
         z$p   <- parnew
         z$ss  <- SSnew
         z$pri <- PPnew
       }
       trypath[[itry+1]] <- z

       if (is.infinite(funpnew))
       { accept <- FALSE
         break
       }
       else {
         alpha <- AlphaFun(trypath[1:(itry+1)])
         dr_steps <- dr_steps +1
         accept <- Accept(alpha)
       }
     }
  }
#-------------------------------------------------------------------------------
# saving if accepted or during burnin in case of adaptive metropolis...
#-------------------------------------------------------------------------------
   if (accept) {
     naccepted <- naccepted +1
     parold    <- parnew           # new parameter set replaces the old one...
     funpold   <- funpnew
   }

# -----  during burnin: SaveIni will be true
   if (SaveIni) {  # use burnin to update covariance...
       ipos <- ipos+1
       if (ipos > outputlength) {
          ipos <- 1
          icov <- outputlength
       }
       pars[,ipos]  <-parold

       if (accept) naccsave2 <- naccsave2+1

       if (i> burnupdate & naccsave2>5 ){  #  5 is arbitrary...
         jj <- max(ipos,icov)
         if (npar > 1) Covar <- (cov(t(pars[,1:jj]))  + diag(1e-16,npar)) * covscale
         else Covar <- var(pars[,1:jj])  + 1e-16 * covscale
         RR    <- try(chol(Covar))
         if (is.numeric(RR)){
           R <- RR
           NewPars    <- NewParsMN   # this function becomes the default update function
           burnupdate <- i + updatecov
           Rnew  <- Rnew + 1
         }
       }
    }
# -----  past burnin  ------

  # current parameter set saved every nupdate steps ...
  if (i == ou1)  {
    SaveIni    <- FALSE
    ii         <-ii+1
    pars[,ii]  <-parold
    funpars[ii]<-funpold
    if (accept) naccsave<- naccsave+1
    ou1 <- i + ou

    if (useSigma)
      sig[ii,]<- 1/sigold

    # update jump parameters by estimating param covariance.. but only if > 5 accepted ones!
    if (isR & ii > nupdate & naccsave > 5) {
       if (npar > 1) Covar <- (cov(t(pars[,1:ii]))  + diag(1e-16,npar)) * covscale
       else Covar <- var(pars[,1:ii])  + 1e-16 * covscale
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

  count <- c(dr_steps=dr_steps, Alfasteps = A_count,num_accepted=naccepted,
             num_covupdate=Rnew)

  res <-list(pars=t(pars),funp=funpars,naccepted=naccepted,sig=sig,
             bestpar= bestPar,bestfunp=bestfunp,count=count)

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
Res <-
  data.frame(rbind(
  mean=apply(object$pars,2,mean),
  sd  =apply(object$pars,2,sd),
  min =apply(object$pars,2,min),
  max =apply(object$pars,2,max),
  q025 =apply(object$pars,2,quantile,probs=0.25),
  q050 =apply(object$pars,2,quantile,probs=0.5),
  q075 =apply(object$pars,2,quantile,probs=0.75)
  ))

if (!is.null(object$sig))
Res <- data.frame(Res,
 sig=rbind(
  mean=apply(object$sig,2,mean),
  sd  =apply(object$sig,2,sd),
  min =apply(object$sig,2,min),
  max =apply(object$sig,2,max),
  q025 =apply(object$sig,2,quantile,probs=0.25),
  q050 =apply(object$sig,2,quantile,probs=0.5),
  q075 =apply(object$sig,2,quantile,probs=0.75)
  ))

return(Res)
}

