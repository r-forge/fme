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
grvar <- expand.grid(map,sensvar)
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
} else {
# global summaries
L1   <- colMeans(abs(Sens))
L2   <- sqrt(colSums(Sens*Sens))/nout
Mean <- colMeans(Sens)
Min  <- apply(Sens,2,min)
Max  <- apply(Sens,2,max)
out  <- data.frame(cbind(value=pp,scale=parscale,L1,L2,Mean,Min,Max))
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

plot.sensFun<- function(x,main=NULL,legpos="topleft",...)
{
  nx  <-attr(x,"nx")
  var <-attr(x,"var")
  Main <- main

 for (i in 1:length(var))
 {
 ii <- ((i-1)*nx):(i*nx)
 if (is.null(main)) Main <- var[i]
 sens<- x[ii,]
 matplot(sens$x,as.matrix( sens[,-(1:2)]),type="l",ylab="-",
        main=Main,...)
 }
  nc <- ncol(x) - 2
  if (! is.na(legpos)) legend(legpos,names(x[,-(1:2)]),col=1:nc,lty=1:nc)

}
