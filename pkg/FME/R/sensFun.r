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

sensFun <- function(parms, # model parameters
                    sensvar=NULL,        # sensitivity variables, default = ALL output variables
                    senspar=names(parms), # sensitivity parameters, default = ALL parameters
                    varscale = NULL, # weighing factor of sensitivity variables, NA = variable value
                    parscale = NULL, # weighing factor of sensitivity parameters, NA = parameter value
                    tiny=1e-8,     # numerical difference factor, perturbation factor
                    solver="ode",  # numerical solution method, one of "ode" or "steady"
                    ...)
{
#----------------------------
# 1. The solver
#----------------------------
if (is.function(solver))
  Solve <- function(parms) solver(parms,...) else
if (solver == "ode")
  Solve <- function(parms) ode(parms=parms,...) else
if (solver == "ode.band")
  Solve <- function(parms) ode.band(parms=parms,...) else
if (solver == "ode.1D")
  Solve <- function(parms) ode.1D(parms=parms,...) else
if (solver == "ode.2D")
  Solve <- function(parms) ode.2D(parms=parms,...) else
if (solver == "steady")   # 1 row, keep variable names...
  Solve <- function(parms)
   { res<-unlist(steady(parms=parms,...))
    cn <- names(res)
    matrix(nr=1,data=res,dimnames=list(NULL,cn))}      else
if (solver == "steady.1D")
  Solve <- function(parms)
   {res <- unlist(steady.1D(parms=parms,...))
    cn <- names(res)
    matrix(nr=1,data=res,dimnames=list(NULL,cn))} else
if (solver == "steady.band")
  Solve <- function(parms)
   { res <- unlist(steady.band(parms=parms,...))
    cn <- names(res)
    matrix(nr=1,data=res,dimnames=list(NULL,cn))} else
if (solver == "steady.2D")
  Solve <- function(parms)
  { res<- unlist(steady.2D(parms=parms,...))
    cn <- names(res)
    matrix(nr=1,data=res,dimnames=list(NULL,cn))} else
stop("Cannot proceed: solver not known ")

yRef  <- Solve(parms)

# if a data.frame or a vector is returned, make it a matrix
if (is.data.frame(yRef)) yRef <- as.matrix(yRef)
if (is.vector(yRef)) yRef<- matrix(nr=1,yRef)

map   <- 1:nrow(yRef)    # e.g. time, the descriptor (x)-variable

#----------------------------
# 2. sensitivity variables
#----------------------------
if (is.null(sensvar))
 {
  ivar    <- 1:ncol(yRef)
  if (is.character(solver))
    if (solver %in% c("ode","ode.band","ode.1D","ode.2D"))
    { ivar <- ivar[-1]
      map  <- yRef[,1]
    }
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

nout  <- length(ivar)
ndim  <- nrow(yRef)
grvar <- expand.grid(map,sensvar)
if (ndim ==1) svar <- sensvar else svar <- paste(grvar[,2],grvar[,1],sep="")

yRef  <- as.vector(yRef[,ivar])

#----------------------------
# 3. sensitivity parameters/initial conditions
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

# global summaries
L1   <- colMeans(abs(Sens))
L2   <- sqrt(colSums(Sens*Sens))/nout
Mean <- colMeans(Sens)
Min  <- apply(Sens,2,min)
Max  <- apply(Sens,2,max)
Sns      <- data.frame(cbind(value=pp,scale=parscale,L1,L2,Mean,Min,Max))
rownames(Sns) <- names(pp)

Fun <- NULL
summvar<-NULL
colnames(Sens) <- names(pp)
Sens <- data.frame(x=grvar[,1],var=as.character(grvar[,2]),Sens)
    
# summaries per variable, only if there are more instances of one variable
if (ndim > 1){
summvar <- data.frame(
           L1=unlist(aggregate(abs(Sens[,-(1:2)]),by=list(Sens[,2]),FUN=mean)[,-1]),
           L2=unlist(aggregate(Sens[,-(1:2)]*Sens[,-(1:2)],by=list(Sens[,2]),FUN=sum)[,-1]),
           Mean=unlist(aggregate(Sens[,-(1:2)],by=list(Sens[,2]),FUN=mean)[,-1]),
           Min=unlist(aggregate(Sens[,-(1:2)],by=list(Sens[,2]),FUN=min)[,-1]),
           Max=unlist(aggregate(Sens[,-(1:2)],by=list(Sens[,2]),FUN=max)[,-1]),
           N=unlist(aggregate(Sens[,-(1:2)],by=list(Sens[,2]),FUN=length)[,-1])
           )
summvar$L2 <- sqrt(summvar$L2)/summvar$N
#summvar <- cbind(var=unlist(aggregate(Sens[,-(1:2)],by=list(Sens[,2]),FUN=length)[,1]),summvar)
}
return(list(model=Sns,var=summvar,fun=Sens))
}

