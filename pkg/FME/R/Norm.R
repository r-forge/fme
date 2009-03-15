## -----------------------------------------------------------------------------
## MultiNormal distribution
## -----------------------------------------------------------------------------

Norm <- function(parMean, parCovar, parRange=NULL, num) {

  nc  <- NCOL(parCovar)
  if (nc != length(parMean))
    stop ("cannot generate Normal distribution: parCovar and parMean not compatible")

  if (nc ==1)
    return (matrix(nc=nc,rnorm(n=num,mean=parMean,sd=sqrt(parCovar))))
  R   <- chol(parCovar)   # Cholesky decomposition
  Z   <- matrix(nc=nc, rnorm(num*nc))

  parset<- Z%*%R
  for (i in 1:nc) parset[,i] <- parset[,i]+parMean[i]
    if (! is.null (parRange)) {
      for (i in 1: ncol(parset)) {
        parset[,i] <- pmin(parset[,i],parRange[i,2])
        parset[,i] <- pmax(parset[,i],parRange[i,1])
      }
    }
  colnames(parset) <- names(parMean)
  return(parset)
}