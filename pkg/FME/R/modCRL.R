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
hist.modCRL<- function(x, what = 1:ncol(x),...) {
  hh <- list()
  hh$pars <- x
  hist.modMCMC(hh, what=what, Full=FALSE, ...)
}

## -----------------------------------------------------------------------------
pairs.modCRL<- function(x, what = 1:ncol(x), nsample = NULL, ...) {

  panel.cor <- function(x, y,...)
    text(x = mean(range(x)), y = mean(range(y)),
        labels = format(cor(x, y), digits = 2))
  panel.hist <- function(x,...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 2))
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "grey")
  }
  panel.main <- function(x,y,...)
    points(x[ii],y[ii],...)

  var <- colnames(x)
  if (! is.numeric(what)) {
      ln <- length(what)
      Select <- which (var %in% what)
      if(length(Select) != ln)
        stop("not all parameters in 'what' are in 'x'")
      what <- Select
  } else {
      if (max(what) > ncol(x))
        stop("index in 'what' too large")
      if (min(what) < 1)
        stop("index in 'what' should be >0")
  }
  X <- x[,what]

  X<- as.matrix(X)

  if (is.null(nsample))
    ii <- 1:nrow(X) else
    ii <- sample((1:nrow(X)),nsample)

  labels <- colnames(X)
  dots <- list(...)

  dots$diag.panel <- if(is.null(dots$diag.panel)) panel.hist else dots$diag.panel
  dots$lower.panel <- if(is.null(dots$lower.panel)) panel.cor else dots$lower.panel
  dots$upper.panel <- if(is.null(dots$upper.panel)) panel.main else dots$upper.panel
  dots$gap <- if(is.null(dots$gap)) 0 else dots$gap
  dots$labels <- if(is.null(dots$labels)) labels else dots$labels

  do.call("pairs",c(alist(X),dots))
}

## -----------------------------------------------------------------------------
plot.modCRL<-function(x, what=NULL, trace = TRUE, ...) {
  vars <- attr(x,"var")

  var <- colnames(x)
  NP  <- attr(x,"npar")
  parnames <- var[1:NP]
  varnames <- var[-(1:NP)]
  if (!is.null(what)) {
    if (! is.numeric(what)) {
      ln <- length(what)
      Select <- NP+which (varnames %in% what)
      SelPar <- which (parnames %in% what)
      if(length(Select) +length(SelPar)!= ln)
        stop("not all variables or parameters in 'what' are in 'x'")
    } else {
      Select <- what[what> NP]
      SelPar <- what[what<= NP]
      if (max(what) > length(var))
        stop("index in 'what' too large")
    }
  } else {
    Select <- NP+(1:length(varnames))
    SelPar <- 1:NP
  }
  nvar <- length(Select)
  np   <- length(SelPar)
  if (np == 0)
    stop("cannot plot monte carlo results: select at least one parameter")
  if (nvar == 0)
    stop("cannot plot monte carlo results: select at least one variable")

  if (np > 1 & nvar >1)
    stop("cannot plot monte carlo results for more than 1 parameter and more than one variable")
  pars <- attr(x,"names")

  dots <- list(...)
  nmdots <- names(dots)

  nv <- max(nvar,np)
  if (! "mfrow" %in% nmdots) {

    nc <- ceiling(sqrt(nv))
    nr <- ceiling(nv/nc)
    mfrow <- c(nr,nc)
  } else mfrow <- dots$mfrow

  if (! is.null(mfrow)) {
    mf <- par(mfrow=mfrow)
#    on.exit(par(mf))
  }
  Main <- is.null(dots$main)
  Xlab <- is.null(dots$xlab)
  dots$ylab    <- if(is.null(dots$ylab))  ""     else dots$ylab

  if (np == 1) {
    ip <- SelPar[1]
    dots$xlab    <- if(is.null(dots$xlab))  pars[ip]   else dots$xlab

    for(i in Select) {
      if (Main) dots$main <- var[i]
      do.call("plot",c(alist(x[,ip],x[,i]),dots))
      if (trace) lines(lowess(x[,ip],x[,i]),lwd=2)
    }
  } else {  # nvar==1
    iv <- Select[1]
    for(i in SelPar) {
      if (Xlab) dots$xlab <- pars[i]
      if (Main) dots$main <- var[iv]
      do.call("plot",c(alist(x[,i],x[,iv]),dots))
      if (trace) lines(lowess(x[,i],x[,iv]),lwd=2)
    }

  }
}
