
###############################
# Model sensitivity
###############################

sensRange <- function(  func,
                        parms, # model parameters
                        sensvar=NULL,   # sensitivity variables, default = ALL output variables
                        dist="unif",    # type of distribution, one of "norm", "unif", "latin", "grid"
                        parRange=NULL,  # if unif: parameter ranges (min, max)
                        parMean=NULL,   # if normal / multinormal: mean value of parameters
                        parCovar=NULL,  # if normal / multinormal: parameter covariances
                        map = 1,        # column number with "mapping variable"
                        num = 100,      # number of runs (approximate if dist="grid"...
                        ...)
{
  Solve <- function(parms) func(parms,...)

if (is.matrix(parms))
  {
  dist<-"input"
  nr <- nrow(parms)
  if (num >= nr) ii <- 1: num else
                 ii <- sample(1:nr,size=num,replace=FALSE)
  parset <-parms[ii,]
  Parms <- parset[1,]
   } else Parms <- parms

# reference run
yRef  <- Solve(Parms)

if (is.vector(yRef)) yRef<- matrix(nr=1,yRef)
if (is.data.frame(yRef)) yRef <- as.matrix(yRef)

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

yRef  <- as.vector(yRef[,ivar])
Sens  <- matrix(ncol=length(yRef),nrow=num,NA)

# sensitivity parameters
senspar <- NULL
senspar <- names(parMean)
if (is.null(senspar)) senspar <- rownames(parRange)
if (is.null(senspar)) senspar <- rownames(parCovar)
if (is.null(senspar)) senspar <- colnames(parCovar)
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
  Parms[ipar]  <- parset[i,]
  yRef <- Solve(Parms)
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
x  <- attr(object,"x")
SumSens <- data.frame(
x    = x,
Mean = apply(sens,2,FUN=mean),
Sd   = apply(sens,2,FUN=sd),
Min  = apply(sens,2,FUN=min),
Max  = apply(sens,2,FUN=max),
q25  = apply(sens,2,FUN=function(x)quantile(x,probs=0.25)),
q50  = apply(sens,2,FUN=function(x)quantile(x,probs=0.5)),
q75  = apply(sens,2,FUN=function(x)quantile(x,probs=0.75))
)
attr(SumSens,"var") <- attr(object,"var")
attr(SumSens,"nx")  <- attr(object,"nx")
class(SumSens)<-c("summary.sensRange","data.frame")

return(SumSens)
}

plot.sensRange<-function(x,main=NULL,...)
{
npar <- attr(x,"npar")
nx  <-attr(x,"nx")
var <-attr(x,"var")
Main <- main

for (i in 1:length(var)){
  if (is.null(main)) Main <- var[i]
  ii <- ((i-1)*nx):(i*nx)
  sens<- x[ii,-(1:npar)]
  matplot(x[ii,2],t(sens),type="l",main=Main,...)
}
}

plot.summary.sensRange<-function(x,main=NULL,legpos="topleft",...)
{
nx  <-attr(x,"nx")
var <-attr(x,"var")
Main <- main
for (i in 1:length(var)){
  ii <- ((i-1)*nx+1):(i*nx)
  X<- x[ii,]
  yrange<-(range(cbind(X$Min,X$Max)))
  if (is.null(main)) Main <- var[i]
  plot(X$x,X$Mean,ylim=yrange,type="n",main=Main,...)
  polygon(c(X$x,rev(X$x)),c(X$Min,rev(X$Max)),col=grey(0.8),border=NA)
  polygon(c(X$x,rev(X$x)),c(X$Mean-X$Sd,rev(X$Mean+X$Sd)),col=grey(0.7),border=NA)
  lines(X$x,X$Mean,lwd=2)
 }
if (! is.null(legpos))
legend(legpos,fill=c(grey(0.9),grey(0.8)),
       legend=c("Min-Max","Mean+-sd"),bty="n")

}
