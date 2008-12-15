##########################################################
Norm <- function(parMean,       # the mean for each parameter
                 parCovar,      # parameter variance-covariance matrix
                 parRange=NULL, # the min and max for each parameter
                 num)           # the number of parameter sets

#--------------------------------------------------------
# draws a multinormally distributed sample; input is a
# covariance matrix and mean
#---------------------------------------------------------

{
nc  <- NCOL(parCovar)
if (nc != length(parMean))
 stop ("cannot generate Normal distribution: parCovar and parMean not compatible")

if (nc ==1) return (matrix(nc=nc,rnorm(n=num,mean=parMean,sd=sqrt(parCovar))))
R   <- chol(parCovar)   # Cholesky decomposition
Z   <- matrix(nc=nc, rnorm(num*nc))

parset<- Z%*%R
for (i in 1:nc) parset[,i] <- parset[,i]+parMean[i]
  if (! is.null (parRange)) {
   for (i in 1: ncol(parset))
    {parset[,i] <- pmin(parset[,i],parRange[i,2])
     parset[,i] <- pmax(parset[,i],parRange[i,1]) }
  }
colnames(parset) <- names(parMean)
return(parset)
}


##########################################################

Latinhyper <- function(parRange,  # the min and max for each parameter
                       num)       # the number of parameter sets
#---------------------------------------------------------------------------
# Draws a latin hypercube sample from parameters whose ranges are given.
#---------------------------------------------------------------------------

{

npar <- nrow(parRange)
latin <- NULL
for (i in 1:npar)
{
 pr    <- unlist(parRange[i,])                     # min and max of this parameter
 dpar  <- diff(pr)/num                             # delt of parameter interval
 ii    <- sort(runif(num),index.return=TRUE)$ix-1  # index to interval (0=1st interval)
 rr    <- runif(num)                               # random number per ii to take a value within each interval
 pval  <- pr[1]+ii*dpar+rr*dpar                    # parameter values
 latin <- cbind(latin,pval)
}
colnames(latin) <- rownames(parRange)
latin
}

##########################################################

Grid <- function(parRange,  # the min and max for each parameter
                 num)       # the number of parameter sets
#---------------------------------------------------------------------------
# Generates a parameter set arranged on a grid
#---------------------------------------------------------------------------

{
npar <- nrow(parRange)
ngrid <- trunc(num^(1/npar))   # number of grid cells per parameter
pargrid <- list()

for (i in 1:npar)
{
 pr   <- unlist(parRange[i,])                  # min and max of this parameter
 pval <- seq(from=pr[1],to=pr[2],len=ngrid)   # parameter values
 pargrid[[i]] <- pval
}
pg <- as.matrix(expand.grid(pargrid))
colnames(pg) <- rownames(parRange)

return(pg)
}

##########################################################

Unif <- function(parRange,  # the min and max for each parameter
                 num)       # the number of parameter sets
#---------------------------------------------------------------------------
# Draws a uniform random sample from a set of parameters
#---------------------------------------------------------------------------

{

npar <- nrow(parRange)

parset   <- NULL
     for (i in 1:npar) parset <- cbind(parset,
                        runif(num,parRange[i,1],parRange[i,2]))
colnames(parset) <- rownames(parRange)
parset
}

