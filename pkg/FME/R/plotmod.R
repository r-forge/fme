plotmod <- function (mat,x=1,which=2:ncol(mat), trace=FALSE, ...) {

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

  if (! is.numeric(which)) {
      ln <- length(which)
      Select <- which (var %in% which)
      if(length(Select) != ln)
        stop("not all parameters in 'which' are in 'x$pars'")
      which <- Select
  } else {
      if (max(which) > ncol(mat))
        stop("index in 'which' too large")
      if (min(which) < 1)
        stop("index in 'which' should be >0")
  }

  np <- length(which)

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

  for(i in which) {
    if (Main) dots$main <- colnames(mat)[i]
    do.call("plot",c(alist(mat[,x],mat[,i]),dots))
    if (trace) lines(lowess(mat[,i]),col="darkgrey",lwd=2)
  }
}

histmod <- function (mat,which=2:ncol(mat),  ...) {

  var <- colnames(mat)
  if (! is.numeric(which)) {
      ln <- length(which)
      Select <- which (var %in% which)
      if(length(Select) != ln)
        stop("not all parameters in 'which' are in 'x$pars'")
      which <- Select
  } else {
      if (max(which) > ncol(mat))
        stop("index in 'which' too large")
      if (min(which) < 1)
        stop("index in 'which' should be >0")
  }

  np <- length(which)

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

    for(i in which) {
    if (Main) dots$main <- colnames(mat)[i]
    do.call("hist",c(alist(mat[,i]),dots))
  }
}
