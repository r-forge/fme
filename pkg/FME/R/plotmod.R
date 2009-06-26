plotmod <- function (mat,x=1,what=2:ncol(mat), trace=FALSE, ...) {

  var <- colnames(mat)
  if (length(x) != 1) stop ("x should contain one value")
  
  if (! is.numeric(x)) {
      Select <- which (var %in% x)
      x <- Select
  } else {
      if (max(x) > ncol(mat))
        stop("index in 'x' too large")
      if (min(x) < 1)
        stop("index in 'x' should be >0")
  }

  if (! is.numeric(what)) {
      ln <- length(what)
      Select <- which (var %in% what)
      if(length(Select) != ln)
        stop("not all parameters in 'what' are in 'x$pars'")
      what <- Select
  } else {
      if (max(what) > ncol(mat))
        stop("index in 'what' too large")
      if (min(what) < 1)
        stop("index in 'what' should be >0")
  }

  np <- length(what)

  dots <- list(...)
  nmdots <- names(dots)

  if (! "mfrow" %in% nmdots) {
    nc <- ceiling(sqrt(np))
    nr <- ceiling(np/nc)
    mfrow <- c(nr,nc)
  } else mfrow <- dots$mfrow

  if (! is.null(mfrow)) {
    mf <- par(mfrow=mfrow)
#    on.exit(par(mf))
  }

  Main <- is.null(dots$main)

  dots$xlab <- if(is.null(dots$xlab)) colnames(mat)[x]else dots$xlab
  dots$ylab <- if(is.null(dots$ylab)) ""  else dots$ylab

  for(i in what) {
    if (Main) dots$main <- colnames(mat)[i]
    do.call("plot",c(alist(mat[,x],mat[,i]),dots))
    if (trace) lines(lowess(mat[,i]),col="darkgrey",lwd=2)
  }
}

histmod <- function (mat,what=2:ncol(mat),  ...) {

  var <- colnames(mat)
  if (! is.numeric(what)) {
      ln <- length(what)
      Select <- which (var %in% what)
      if(length(Select) != ln)
        stop("not all parameters in 'what' are in 'x$pars'")
      what <- Select
  } else {
      if (max(what) > ncol(mat))
        stop("index in 'what' too large")
      if (min(what) < 1)
        stop("index in 'what' should be >0")
  }

  np <- length(what)

  dots <- list(...)
  nmdots <- names(dots)

  if (! "mfrow" %in% nmdots) {
    nc <- ceiling(sqrt(np))
    nr <- ceiling(np/nc)
    mfrow <- c(nr,nc)
  } else mfrow <- dots$mfrow

  if (! is.null(mfrow)) {
    mf <- par(mfrow=mfrow)
#    on.exit(par(mf))
  }

  Main <- is.null(dots$main)

  dots$ylab <- if(is.null(dots$ylab)) ""  else dots$ylab

    for(i in what) {
    if (Main) dots$main <- colnames(mat)[i]
    do.call("hist",c(alist(mat[,i]),dots))
  }
}
