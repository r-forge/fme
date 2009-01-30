
###############################
# Model sensitivity
###############################

sensRange <- function(  func,
                        parms=NULL, # model parameters
                        sensvar=NULL,   # sensitivity variables, default = ALL output variables
                        dist="unif",    # type of distribution, one of "norm", "unif", "latin", "grid"
                        parInput=NULL,  # random parameter values
                        parRange=NULL,  # if unif: parameter ranges (min, max)
                        parMean=NULL,   # if normal / multinormal: mean value of parameters
                        parCovar=NULL,  # if normal / multinormal: parameter covariances
                        map = 1,        # column number with "mapping variable"
                        num = 100,      # number of runs (approximate if dist="grid"...
                        ...)
{
if (is.null(parms) & ! is.null(parInput))
  parms <- parInput[1,]
if (is.null(parms)) parms <- parMean
if (is.null(parms)) stop ("'parms' not known")
if (is.matrix(parms) && nrow(parms)>1) stop ("'parms' should be a vector")

Solve <- function(parms) func(parms,...)

if (! is.null(parInput))
  {
  dist<-"input"
  nr <- nrow(parInput)
  num <- min(num,nr)
  if (num == nr) ii <- 1: nr else
                 ii <- sample(1:nr,size=num,replace=FALSE)
  parset <-as.matrix(parInput[ii,])
  if(is.null(parms)) parms <- parInput[1,]
  }

 Parms <- parms

# reference run
yRef  <- Solve(Parms)

if (is.vector(yRef))
  {
    ynames <- names(yRef)
    yRef <- matrix(nr=1,yRef)
    colnames(yRef) <- ynames
  } else if (is.data.frame(yRef)) yRef <- as.matrix(yRef)

# check sensitivity variables
if (is.null(sensvar))
 {
  ivar    <- 1:ncol(yRef)
    if (! is.null(map)) ivar <- ivar[-map]
  sensvar <- colnames(yRef)[ivar]
  if(is.null(sensvar))  sensvar <- ivar
 } else {
  ivar  <- findvar(yRef[1,],sensvar,"variables")
  if (! is.character(sensvar)) # try to create names rather than nrs
  {sv <- sensvar
    sensvar<-colnames(yRef)[ivar]
   if (is.null(sensvar)) sensvar<-sv
  }
 }
if (is.null(map)) map   <- 1:nrow(yRef) else map <- yRef[,map]

nout  <- length(ivar)
ndim  <- nrow(yRef)
grvar <- expand.grid(map,sensvar)
if (ndim ==1) svar <- sensvar else svar <- paste(grvar[,2],grvar[,1],sep="")

YREF  <- as.vector(yRef[,ivar])
Sens  <- matrix(ncol=length(YREF),nrow=num,NA)

# sensitivity parameters
senspar <- NULL
senspar <- names(parMean)
if (is.null(senspar)) senspar <- rownames(parRange)
if (is.null(senspar)) senspar <- rownames(parCovar)
if (is.null(senspar)) senspar <- colnames(parCovar)
if (is.null(senspar)) senspar <- colnames(parInput)
if (is.null(senspar)) senspar <- names(Parms)
#if (is.null(senspar)) stop("parameter names are not known")

npar  <- length(senspar)

ipar <- findvar(Parms,senspar,"parameters")
pp   <- unlist(Parms)[ipar]

# sanity checks for random parameters
if (dist == "norm" && (is.null(parMean) | is.null(parCovar)))
  stop("parMean and parCovar should be given if dist = norm")
if(!(dist  %in% c("norm","input")) && is.null(parRange))
  stop("parRange should be given if dist = unif, grid or latin")

# generate random parameters
if (dist == "norm")
 parset <- Norm (parCovar=parCovar, parMean=parMean, parRange, num) else
if (dist == "unif")
 parset <- Unif(parRange,num)                      else
if (dist =="latin")
 parset <- Latinhyper(parRange,num)                else
if (dist =="grid")
 parset <- Grid(parRange,num)                      else
if (dist != "input" )stop("dist should be one of 'norm','unif','latin' or 'grid' ")


# The sensitivity output
colnames(Sens) <- svar

for (i in 1:num)
{
  if(prod(Parms[ipar] == parset[i,])==0) { # no need to run model again if same parameter value (e.g. MCMCrun!)
    Parms[ipar]  <- parset[i,]
    yRef <- Solve(Parms)
  }
  if (is.vector(yRef)) Sens[i,] <- yRef[ivar] else
                       Sens[i,] <- as.vector(unlist(yRef[,ivar]))   # unlist in case it is a data.frame
}
sens=data.frame(cbind(parset,Sens))
class(sens) <- c("sensRange","data.frame")
attr(sens,"npar") <- ncol(parset)
attr(sens,"x")   <- map
attr(sens,"nx")  <- length(map)
attr(sens,"var") <- sensvar
return (sens)
}

summary.sensRange<-function(object,...)
{
npar<-attr(object,"npar")
sens<- as.matrix(object[,-(1:npar)])
x   <- attr(object,"x")
names(x) <- NULL
nx  <-attr(object,"nx")
var <-attr(object,"var")

SumSens <- data.frame(
x    = x,
Mean = apply(sens,2,FUN=mean),
Sd   = apply(sens,2,FUN=sd),
Min  = apply(sens,2,FUN=min),
Max  = apply(sens,2,FUN=max),
q05  = apply(sens,2,FUN=function(x)quantile(x,probs=0.05)),
q25  = apply(sens,2,FUN=function(x)quantile(x,probs=0.25)),
q50  = apply(sens,2,FUN=function(x)quantile(x,probs=0.5)),
q75  = apply(sens,2,FUN=function(x)quantile(x,probs=0.75)),
q95  = apply(sens,2,FUN=function(x)quantile(x,probs=0.95))
)

rownames(SumSens) <- colnames(sens)
attr(SumSens,"var") <- attr(object,"var")
attr(SumSens,"nx")  <- attr(object,"nx")
class(SumSens)<-c("summary.sensRange","data.frame")

return(SumSens)
}


plot.sensRange<-function(x,xyswap=FALSE,what=NULL,...)
{
npar <- attr(x,"npar")
nx  <-attr(x,"nx")
var <-attr(x,"var")
X  <-attr(x,"x")
dots <- list(...)
nmdots <- names(dots)

sens<- x[,-(1:npar)]

if (!is.null(what)) Select <- which (var %in% what) else Select <- 1:length(var)

if (nx > 1)
  for (i in Select){
    yrange<- if ("ylim" %in% nmdots) dots$ylim else range(X)  # USED for swap!
    ii <- ((i-1)*nx+1):(i*nx)

    if ("main" %in% nmdots) {
    if (!xyswap) matplot(X,t(sens[,ii]),type="l",...)
    else matplot(t(sens[,ii]),X,type="l",ylim=rev(yrange))
    }else {
    Main  <- var[i]
    if (!xyswap) matplot(X,t(sens[,ii]),type="l",main=Main,...)
    else matplot(t(sens[,ii]),X,type="l",main=Main,ylim=rev(yrange))
    }
  }
else
 boxplot(sens[,Select])
}

plot.summary.sensRange<-function(x,xyswap=FALSE,
                                 what=NULL,legpos="topleft",
                                 col=c(grey(0.8),grey(0.7)),
                                 quant=FALSE,
                                 ...)
{
nx  <-attr(x,"nx")
var <-attr(x,"var")

dots <- list(...)
nmdots <- names(dots)

if (!is.null(what)) Select <- which (var %in% what) else Select <- 1:length(var)

if (nx > 1)  {    # summary of a times series or a profile...
 for (i in Select){
  ii <- ((i-1)*nx+1):(i*nx)
  X<- x[ii,]
  if (quant)
  {
  xmin <- X$q05
  xmax <- X$q95
  xmean <- X$q50
  xlow <- X$q25
  xup  <- X$q75
  leg  <- c("q05-q95","q25-q75")
  }else
  {
  xmin <- X$Min
  xmax <- X$Max
  xmean <- X$Mean
  xlow <- X$Mean-X$Sd
  xup  <- X$Mean+X$Sd
  leg <- c("Min-Max","Mean+-sd")
  }
  yrange<- if ("ylim" %in% nmdots) dots$ylim else range(cbind(xmin,xmax))
  xrange<- if ("xlim" %in% nmdots) dots$xlim else range(X$x)
  Main  <- if ("main" %in% nmdots) dots$main else var[i]
  xlab  <- if ("xlab" %in% nmdots) dots$xlab else "x"
  ylab  <- if ("ylab" %in% nmdots) dots$ylab else "y"
  if (!xyswap) {
    plot(X$x,xmean,ylim=yrange,xlim=xrange,type="n",main=Main,xlab=xlab,ylab=ylab)
    polygon(c(X$x,rev(X$x)),c(xmin,rev(xmax)),col=col[1],border=NA)
    polygon(c(X$x,rev(X$x)),c(xlow,rev(xup)),col=col[2],border=NA)
    lines(X$x,xmean,...) }
  else {
    plot(xmean,X$x,xlim=yrange,type="n",main=Main,ylim=rev(xrange),xlab=xlab,ylab=ylab)
    polygon(c(xmin,rev(xmax)),c(X$x,rev(X$x)),col=col[1],border=NA)
    polygon(c(xlow,rev(xup)),c(X$x,rev(X$x)),col=col[2],border=NA)
    lines(xmean,X$x,...)
  }
 }
if (! is.null(legpos))
legend(legpos,fill=c(grey(0.9),grey(0.8)),
       legend=leg,bty="n")
} else             # one summary per variable
{
 X <- x[Select,]
 dotchart(X$Mean,labels=rownames(X),xlim=range(c(X$Min,X$Max)),...)
# add ranges
 nr <- nrow(X)

for (i in 1:nr)
 {
 segments(X$Min[i],i,X$Max[i],i,lty=1)
 segments(X$q05[i],i,X$q95[i],i,lty=1)
 segments(X$q25[i],i,X$q75[i],i,lwd=3)
 }

 if (! is.null(legpos))
 legend(legpos,lwd=1:3,legend=c("Min-Max","q25-q75","q05-q95"),bty="n")

}

}