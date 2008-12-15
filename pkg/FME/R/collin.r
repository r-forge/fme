
###############################
# Collinearity indices
###############################

collin <- function(sensfun)  # sensitivity functions, as estimated from sensFun(...,)
{
if (colnames(sensfun)[1]=="x" && colnames(sensfun)[2] == "var")
   Sens <- sensfun[,-(1:2)] else Sens<-sensfun

npar <- ncol(Sens)
if (npar > 14)
    warning ("will reduce collinearity estimates: too many combinations")
normSens <- t(t(Sens) / sqrt(colSums(Sens*Sens)))

combin <- function(n, v)  # combinations of n elements from a vector p (length (p) > = n)
{
 if (n == 1)
    matrix(data=v, ncol=1)
 else if (n >= length(v))
    matrix(data=v, nrow=1)
 else rbind(cbind(v[1],combin(n-1,v[-1])),combin(n,v[-1]))
 }

collin <- NULL

pset   <- 1:npar

for (n in 2:npar)
{
  numcomb <- choose(npar,n)
  if (numcomb < 5000)
  {
  cc  <-combin(n,pset)

  for (i in 1:nrow(cc))
  {
   ii    <- cc[i,]
   S     <- normSens[,ii]
   id    <- 1/sqrt(min(eigen(t(S)%*%S)$value))
   psub     <- rep(0,npar)
   psub[ii] <- 1
   collin   <- rbind(collin,c(psub,n,id))
  }
 }
}
colnames(collin)<-c(colnames(Sens),"N","collinearity")
return(as.data.frame(collin ))
}
  
