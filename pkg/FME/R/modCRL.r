## -----------------------------------------------------------------------------
## Monte Carlo runs
## -----------------------------------------------------------------------------

modCRL <- function( func, parms=NULL, sensvar=NULL, dist="unif",
                   parInput=NULL, parRange=NULL, parMean=NULL,
                   parCovar=NULL, num = 100, ...) {

  if (is.null(parms) & ! is.null(parInput))
    parms <- parInput[1,]
  if (is.null(parms))
    parms <- parMean
  if (is.null(parms)) {
    if (is.vector(parRange))
      parms<-mean(parRange)
    else parms <- rowMeans(parRange)
  }
  if (is.null(parms))
    stop ("'parms' not known")
  if (is.matrix(parms) && nrow(parms)>1)
    stop ("'parms' should be a vector")

  Solve <- function(parms) func(parms,...)

  if (! is.null(parInput)) {
    dist<-"input"
    nr <- nrow(parInput)
    num <- min(num,nr)
    if (num == nr)
      ii <- 1: nr
    else ii <- sample(1:nr,size=num,replace=FALSE)

    parset <-as.matrix(parInput[ii,])
    if(is.null(parms)) parms <- parInput[1,]
  }

  Parms <- parms

  # reference run
  yRef  <- func(parms,...)

  if (is.matrix(yRef) | is.data.frame(yRef))
    if(nrow(yRef)>1)
      stop("func should return a vector or a matrix/data.frame with one row")

  sens <- sensRange(func,parms,sensvar,dist,parInput,parRange,parMean,
                  parCovar,map=NULL,num=num,...)

  class(sens) <- c("modCRL","data.frame")
  return (sens)
}

## -----------------------------------------------------------------------------
## S3 methods of modCRL
## -----------------------------------------------------------------------------

summary.modCRL<-function(object,...) {
  SumSens <- summary.sensRange(object,...)
  class(SumSens)<-c("summary.modCRL","data.frame")
  return(SumSens)
}

## -----------------------------------------------------------------------------

plot.modCRL<-function(x,what=NULL,xlab=NULL,ylab=NULL,...) {
  var <-attr(x,"var")
  np  <- attr(x,"npar")
  if (np > 1)
    stop("cannot plot monte carlo results for more than 1 parameter")
  pars <- attr(x,"names")[1]

  dots <- list(...)
  nmdots <- names(dots)

  if (!is.null(what))
    Select <- which (var %in% what)
  else Select <- 1:length(var)

  nv <- length(Select)
  if (! "mfrow" %in% nmdots) {
    nc <- ceiling(sqrt(nv))
    nr <- ceiling(nv/nc)
    mfrow <- c(nr,nc)
  } else mfrow <- dots$mfrow

  if (! is.null(mfrow)) {
    mf <- par(mfrow=mfrow)
    on.exit(par(mf))
  }
  
  if (is.null(xlab))
    xlab<-pars
  if (is.null(ylab))
    ylab<-ylab

  for(i in Select) {
    if ("main" %in% nmdots)
      plot(x[,1],x[,i+1], xlab=xlab,ylab=ylab,...)
    else
      plot(x[,1],x[,i+1],main=var[i],
                      xlab=pars,ylab="",...)
  }
}
