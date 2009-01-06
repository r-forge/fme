
###############################
# Model sensitivity
###############################

sensRange <- function(  parms, # model parameters
                        sensvar=NULL,   # sensitivity variables, default = ALL output variables
                        dist="unif",    # type of distribution, one of "norm", "unif", "latin", "grid"
                        parRange=NULL,  # if unif: parameter ranges (min, max)
                        parMean=NULL,   # if normal / multinormal: mean value of parameters
                        parCovar=NULL,  # if normal / multinormal: parameter covariances
                        solver="ode",   # numerical solution method, one of "ode" or "steady"
                        num = 100,      # number of runs (approximate if dist="grid"...
                        Full=FALSE,
                        ...)
{
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

map   <- 1:nrow(yRef)    # e.g. time, the descriptor (x)-variable

# check sensitivity variables
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

SumSens <- data.frame(
x    = grvar[,1],
var  = as.character(grvar[,2]),
Mean = apply(Sens,2,FUN=mean),
Sd   = apply(Sens,2,FUN=sd),
Min  = apply(Sens,2,FUN=min),
Max  = apply(Sens,2,FUN=max),
q25  = apply(Sens,2,FUN=function(x)quantile(x,probs=0.25)),
q50  = apply(Sens,2,FUN=function(x)quantile(x,probs=0.5)),
q75  = apply(Sens,2,FUN=function(x)quantile(x,probs=0.75))
)

if (Full) return (list(summ=SumSens,sens=cbind(parset,Sens))) else return(list(summ=SumSens))

}

