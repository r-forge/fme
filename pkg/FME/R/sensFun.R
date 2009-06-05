## private function

findvar <- function(var1, var2, str="var") {
  if (is.character(var2[[1]])){
    ivar  <- which (names(var1)%in%var2)
    if (length(ivar)!= length(var2))
      stop(paste("cannot proceed: not all sensitivity", str,"are known"))
    return(ivar)
  } else {
  if (max(var2)>length(var1))
    stop (paste("cannot proceed: index to sensitivity ", str, "too large"))
  return(var2)
  }
}

## -----------------------------------------------------------------------------
## Sensitivity functions
## -----------------------------------------------------------------------------

sensFun <- function(func, parms, sensvar=NULL, senspar=names(parms),
                    varscale = NULL, parscale = NULL,
                    tiny=1e-8, map = 1, ...) {
  ## 1. The solver
  Solve <- function(parms) func(parms,...)

  yRef  <- Solve(parms)
  Type <- 1
  if (class(yRef)=="modCost") {
    Res    <- yRef$residuals
    ynames <- Res$name
    yRef <- cbind(Res$x,Res$mod)
    names(yRef) <- ynames
    Type <- 2
    sensvar <- NULL   # input of sensvar not allowed for modCost type

    Solve <- function(parms) {
      Res<- func(parms,...)$residuals
      cbind(Res$x,Res$mod)
    }

  }
  ## if a data.frame or a vector is returned, make it a matrix
  if (is.data.frame(yRef)) yRef <- as.matrix(yRef)
  if (is.vector(yRef)) {
    ynames <- names(yRef)
    yRef <- matrix(nr=1,yRef)
    colnames(yRef) <- ynames
  }
  
  ## 2. sensitivity variables
  if (is.null(sensvar)) {
    ivar    <- 1:ncol(yRef)
    if (! is.null(map))
      ivar <- ivar[-map]
    sensvar <- colnames(yRef)[ivar]
    if(is.null(sensvar))
      sensvar <- ivar
  } else {
    ivar  <- findvar(yRef[1,],sensvar,"variables")
    if (! is.character(sensvar)) { # try to create names rather than nrs
      sv <- sensvar
      sensvar<-colnames(yRef)[ivar]
      if (is.null(sensvar))
        sensvar<-sv
    }
  }
  if (is.null(map))
    map   <- 1:nrow(yRef)
  else map <- yRef[,map]

  nout  <- length(ivar)
  ndim  <- nrow(yRef)
  if (Type ==1)
    grvar <- expand.grid(map,sensvar)
  else grvar<-data.frame(x=map,var=ynames)

  if (ndim ==1)
    svar <- sensvar
  else svar <- paste(grvar[,2],grvar[,1],sep="")

  yRef  <- as.vector(yRef[,ivar])

  ## 3. sensitivity parameters/
  npar  <- length(senspar)
  if (npar ==0)
    stop ("cannot proceed: there are no sensitivity parameters")

  ipar <- findvar(parms,senspar,"parameters")
  pp    <- unlist(parms)[ipar]

  ## 4. perturbed parameters
  dp        <- pp*tiny
  dp[dp==0] <- tiny
  ii      <- which (abs(dp)<tiny)
  dp[ii]  <- sign(dp[ii])*tiny

  if (is.null(parscale))
    parscale <- pp
  else parscale<-rep(parscale,npar)

  if (is.null(varscale))
    varscale <- yRef
  else varscale <- rep (varscale,length(yRef))

  ## 0 is set equal to a very small number
  varscale[varscale == 0]<-1e-20
  parscale[parscale == 0]<-1e-20

  Sens    <- matrix(nrow=length(yRef),ncol=npar,NA)

  ## 5. Loop over all parameters
  for (i in 1:length(ipar)) {
    dval    <- pp[i]+dp[i]
    parms[ipar[i]] <- dval
    Yres    <- Solve(parms)
    if (is.vector(Yres))
      yPert <- Yres[ivar]
    else yPert <- as.vector(unlist(Yres[,ivar]))
    Sens[,i]<- (yPert-yRef)/dp[i] *parscale[i] /varscale
    parms[ipar[i]] <- pp[i]
  }

  ## 6. Finally
  colnames(Sens) <- names(pp)
  Sens <- data.frame(x=grvar[,1],var=as.character(grvar[,2]),Sens)
  attr(Sens,"class") <- c("sensFun","data.frame")
  attr(Sens,"pars") <- pp
  attr(Sens,"parscale") <- parscale
  attr(Sens,"varscale") <- varscale
  if (Type==2) {
    attr(Sens,"var") <- as.vector(unique(ynames))
    attr(Sens,"nx") <- c(0,cumsum(as.vector(table(ynames))))  #start of each var
  } else {
    attr(Sens,"var") <- sensvar
    attr(Sens,"nx") <- length(map)
  }
  attr(Sens,"x")   <- map
  attr(Sens,"Type") <- Type  # type 1: modCost
  return(Sens)
}

## -----------------------------------------------------------------------------
## S3 methods of sensFun
## -----------------------------------------------------------------------------

summary.sensFun <- function(object,vars=FALSE,...) {

  pp       <- attributes(object)$pars
  parscale <-  attributes(object)$parscale
  Sens <- object[,-(1:2)]
  nout <- nrow(Sens)
  if (vars) { # summaries per variable
    Vars <- object[,2]
    out <- data.frame(
         L1  =unlist(aggregate(abs(Sens),by=list(Vars),FUN=mean)[,-1]),
         L2  =unlist(aggregate(Sens*Sens,by=list(Vars),FUN=sum)[,-1]),
         Mean=unlist(aggregate(Sens,by=list(Vars),FUN=mean)[,-1]),
         Min =unlist(aggregate(Sens,by=list(Vars),FUN=min)[,-1]),
         Max =unlist(aggregate(Sens,by=list(Vars),FUN=max)[,-1]),
         N   =unlist(aggregate(Sens,by=list(Vars),FUN=length)[,-1])
      )
    out$L2 <- sqrt(out$L2/out$N)
    out$var <- unique(Vars)
    np <- length(pp)
    nv <- length(unique(Vars))
    out <- data.frame(cbind(value=rep(pp,times=rep(nv,np)),
                      scale=rep(parscale,times=rep(nv,np)),out))
  } else {  # global summaries
    L1   <- colMeans(abs(Sens))
    L2   <- sqrt(colSums(Sens*Sens))/nout
    Mean <- colMeans(Sens)
    Min  <- apply(Sens,2,min)
    Max  <- apply(Sens,2,max)
    N  <- apply(Sens,2,length)
    out  <- data.frame(cbind(value=pp,scale=parscale,L1,L2,Mean,Min,Max,N))
    rownames(out) <- names(pp)
  }
  class(out) <- c("summary.sensFun","data.frame")
  return(out)
}

## -----------------------------------------------------------------------------

print.summary.sensFun<-function(x,...)
  print(format(x,digits=2))

## -----------------------------------------------------------------------------

pairs.sensFun <- function (x, what=NULL, ...) {

  if (!is.null(what)) {
    nx  <-attr(x,"nx")
    var <-attr(x,"var")
    TYP <-attr(x,"Type")

    if (! is.numeric(what)) {
      ln <- length(what)
      Select <- which (var %in% what)
      if(length(Select) != ln)
        stop("not all variables in 'what' are in 'x'")
    } else {
       Select <- what
       if (max(Select) > nx)
         stop("index in 'what' too large")
    }
    ii <- NULL
    
    if (TYP == 1)
     for (i in Select)
       ii <- c(ii,((i-1)*nx):(i*nx))
    else  for (i in Select)
       ii <- c(ii,(nx[i]+1):nx[i+1])
    }
  else ii <- 1:nrow(x)

  if (colnames(x)[1]=="x" && colnames(x)[2] == "var")
    X <- x[ii,-(1:2)] else X<-x[ii,]

  panel.cor <- function(x, y) text(x = mean(range(x)), y = mean(range(y)),
       labels = format(cor(x, y), digits = 2))
  dots <- list(...)

  dots$diag.panel <- if(is.null(dots$diag.panel)) NULL else dots$diag.panel
  dots$lower.panel <- if(is.null(dots$lower.panel)) panel.cor else dots$lower.panel
  dots$gap <- if(is.null(dots$gap)) 0 else dots$gap
  do.call("pairs",c(alist(as.matrix(X)),dots))

}

## -----------------------------------------------------------------------------

plot.sensFun<- function(x, what=NULL, legpos="topleft", ...) {
  nx  <-attr(x,"nx")
  var <-attr(x,"var")
  TYP <-attr(x,"Type")

  dots <- list(...)
  nmdots <- names(dots)

  nc <- ncol(x) - 2

  Main <- is.null(dots$main)
  dots$ylab <- if(is.null(dots$ylab)) "sensitivity" else dots$ylab
  dots$type <- if(is.null(dots$type)) "l" else dots$type
  dots$col <- if(is.null(dots$col)) 1:nc else dots$col

  if (!is.null(what)) {
    if (! is.numeric(what)) {

      ln <- length(what)
      Select <- which (var %in% what)
      if(length(Select) != ln)
        stop("not all variables in 'what' are in 'x'")
    } else {     # index
      Select <- what
      if (max(Select) > nx)
        stop("index in 'what' too large")
    }
    if(is.null(dots$ylim)) {
      ii <- NULL
      for (i in Select){
        if (TYP == 1)
          ii <- c(ii,((i-1)*nx):(i*nx))
        else
          ii <- c(ii,(nx[i]+1):nx[i+1])
          dots$ylim <- range(x[ii,-(1:2)])
      }
    }
  }  else {
    Select <- 1:length(var)
    dots$ylim <- if(is.null(dots$ylim)) range(x[,-(1:2)])
  }
  Lty <- is.null(dots$lty)

  st <- 1
  for (i in Select){
    if (TYP == 1)
      ii <- ((i-1)*nx):(i*nx)
    else
      ii <- (nx[i]+1):nx[i+1]
    if (Main) dots$main <- var[i]
    sens<- x[ii,]
    dots$lty <- if(Lty) st
    if (st==1)
      do.call("matplot",c(alist(sens$x,as.matrix( sens[,-(1:2)])),dots))
    else
      do.call("matlines",c(alist(sens$x,as.matrix( sens[,-(1:2)])),dots))
    st <- st+1
  }
  if (! is.na(legpos))
    legend(legpos,names(x[,-(1:2)]),col=1:nc,lty=1)

}

## -----------------------------------------------------------------------------

plot.summary.sensFun<- function(x,...) {
  dots <- list(...)
  nmdots <- names(dots)

  nr <- nrow(x)
  dots$main <- if(is.null(dots$main)) "sensitivity" else dots$main
  dots$labels<- if(is.null(dots$labels)) rownames(x) else dots$labels
  dots$xlim <- if(is.null(dots$xlim)) range(c(x$Min,x$Max)) else dots$xlim
  dots$pch <- if(is.null(dots$pch)) 1 else dots$pch
  dots$col <- if(is.null(dots$col)) "black" else dots$col

  do.call("dotchart",c(alist(c(x$Mean,Inf)),dots))

  # add ranges
  for (i in 1:nr) {
    segments(x$Min[i],i,x$Max[i],i,lty=1)
  }
  points(x$L1,1:nr,pch=16,col="red")
  points(x$L2,1:nr,pch=18,col="blue")

  legend("top",legend=c("L1","L2","Mean"),pch=c(dots$pch ,16,18),
          col=c(dots$col,"red","blue"),ncol=3)
}
