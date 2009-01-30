# private function

findvar <- function(var1,
                    var2,
                    str="var")
{
 if (is.character(var2[[1]])){
  ivar  <- which (names(var1)%in%var2)
  if (length(ivar)!= length(var2)) stop(paste("cannot proceed: not all sensitivity", str,"are known"))
  return(ivar)
 } else
 {if (max(var2)>length(var1)) stop (paste("cannot proceed: index to sensitivity ", str, "too large"))
  return(var2)
  }
}

###############################
# Sensitivity functions
###############################

sensFun <- function(func,
                    parms, # model parameters
                    sensvar=NULL,        # sensitivity variables, default = ALL output variables
                    senspar=names(parms), # sensitivity parameters, default = ALL parameters
                    varscale = NULL, # weighing factor of sensitivity variables, NA = variable value
                    parscale = NULL, # weighing factor of sensitivity parameters, NA = parameter value
                    tiny=1e-8,     # numerical difference factor, perturbation factor
                    map = 1,
                    ...)
{
#----------------------------
# 1. The solver
#----------------------------
Solve <- function(parms) func(parms,...)

yRef  <- Solve(parms)
Type <- 1

if (class(yRef)=="modCost")
  { Res    <- yRef$residuals
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
# if a data.frame or a vector is returned, make it a matrix
if (is.data.frame(yRef)) yRef <- as.matrix(yRef)
if (is.vector(yRef))
  {
    ynames <- names(yRef)
    yRef <- matrix(nr=1,yRef)
    colnames(yRef) <- ynames
  }
  
#----------------------------
# 2. sensitivity variables
#----------------------------
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
if (Type ==1) grvar <- expand.grid(map,sensvar) else
grvar<-data.frame(x=map,var=ynames)
if (ndim ==1) svar <- sensvar else svar <- paste(grvar[,2],grvar[,1],sep="")

yRef  <- as.vector(yRef[,ivar])

#----------------------------
# 3. sensitivity parameters/
#----------------------------
# check sensitivity parameters and stop if not all known
npar  <- length(senspar)
if (npar ==0) stop ("cannot proceed: there are no sensitivity parameters")

ipar <- findvar(parms,senspar,"parameters")
pp    <- unlist(parms)[ipar]


#----------------------------
# perturbed parameters
#----------------------------
dp        <- pp*tiny   #Change if 0!    KARLINE: add initial conditions...
dp[dp==0] <- tiny
ii      <- which (abs(dp)<tiny)
dp[ii]  <- sign(dp[ii])*tiny

if (is.null(parscale)) parscale <- pp else parscale<-rep(parscale,npar)
if (is.null(varscale)) varscale <- yRef else varscale <- rep (varscale,length(yRef))

# 0 is set equal to a very small number (e-10)
varscale[varscale == 0]<-1e-20
parscale[parscale == 0]<-1e-20

Sens    <- matrix(nrow=length(yRef),ncol=npar,NA)

for (i in 1:length(ipar))
{
  dval    <- pp[i]+dp[i]
  parms[ipar[i]] <- dval
  Yres    <- Solve(parms)
  if (is.vector(Yres)) yPert <- Yres[ivar] else
                       yPert <- as.vector(unlist(Yres[,ivar]))
  Sens[,i]<- (yPert-yRef)/dp[i] *parscale[i] /varscale
  parms[ipar[i]] <- pp[i]
}
colnames(Sens) <- names(pp)

Sens <- data.frame(x=grvar[,1],var=as.character(grvar[,2]),Sens)
attr(Sens,"class") <- c("sensFun","data.frame")
attr(Sens,"pars") <- pp
attr(Sens,"parscale") <- parscale
attr(Sens,"varscale") <- varscale
attr(Sens,"var") <- sensvar
attr(Sens,"nx") <- length(map)
attr(Sens,"x")   <- map
attr(Sens,"Type") <- Type  # type 1: modCost
return(Sens)
}


#S3 methods of sensitivity functions


summary.sensFun <- function(object,vars=FALSE,...)
{
pp       <- attributes(object)$pars
parscale <-  attributes(object)$parscale
Sens <- object[,-(1:2)]
nout <- nrow(Sens)
if (vars){
Vars <- object[,2]
out <- data.frame(
           L1  =unlist(aggregate(abs(Sens),by=list(Vars),FUN=mean)[,-1]),
           L2  =unlist(aggregate(Sens*Sens,by=list(Vars),FUN=sum)[,-1]),
           Mean=unlist(aggregate(Sens,by=list(Vars),FUN=mean)[,-1]),
           Min =unlist(aggregate(Sens,by=list(Vars),FUN=min)[,-1]),
           Max =unlist(aggregate(Sens,by=list(Vars),FUN=max)[,-1]),
           N   =unlist(aggregate(Sens,by=list(Vars),FUN=length)[,-1])
           )
out$L2 <- sqrt(out$L2)/out$N
out$var <- unique(Vars)
np<-length(pp)
nv<-length(unique(Vars))
out <- data.frame(cbind(value=rep(pp,times=rep(nv,np)),
                        scale=rep(parscale,times=rep(nv,np)),out))
} else {
# global summaries
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

print.summary.sensFun<-function(x,...)
  {
  print(format(x,digits=2))
  }


pairs.sensFun <- function (x, ...)
{
    if (colnames(x)[1]=="x" && colnames(x)[2] == "var")
    X <- x[,-(1:2)] else X<-x

    panel.cor <- function(x, y) text(x = mean(range(x)), y = mean(range(y)),
        labels = format(cor(x, y), digits = 2))
    pairs(as.matrix(X), diag.panel = NULL, gap = 0,
        lower.panel = panel.cor, ...)
}

plot.sensFun<- function(x,legpos="topleft",...)
{
  nx  <-attr(x,"nx")
  var <-attr(x,"var")
  dots <- list(...)
  nmdots <- names(dots)

  for (i in 1:length(var)){
   ii <- ((i-1)*nx):(i*nx)
   Main  <- if ("main" %in% nmdots) dots$main else var[i]
   sens<- x[ii,]
   matplot(sens$x,as.matrix( sens[,-(1:2)]),type="l",ylab="-",
        main=Main,...)
  }
  nc <- ncol(x) - 2
  if (! is.na(legpos)) legend(legpos,names(x[,-(1:2)]),col=1:nc,lty=1:nc)

}

plot.summary.sensFun<- function(x,...)
{
  dots <- list(...)
  nmdots <- names(dots)

  nr <- nrow(x)

  if ("main" %in% nmdots)
    dotchart(c(x$Mean,Inf),labels=rownames(x),xlim=range(c(x$Min,x$Max)),
             pch=3,...)
  else dotchart(c(x$Mean,Inf),labels=rownames(x),xlim=range(c(x$Min,x$Max)),
             pch=3, main="Sensitivity")

 # add ranges
 for (i in 1:nr)
 {
 segments(x$Min[i],i,x$Max[i],i,lty=1)
 }
 points(x$L1,1:nr,pch=16,col="red")
 points(x$L2,1:nr,pch=18,col="blue")

 legend("top",legend=c("L1","L2","Mean"),pch=c(3,16,18),col=c("black","red","blue"),ncol=3)
}
