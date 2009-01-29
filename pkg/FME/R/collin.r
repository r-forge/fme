
###############################
# Collinearity indices
###############################

collin <- function(sensfun,parset=NULL)  # sensitivity functions, as estimated from sensFun(...,)
{
if (colnames(sensfun)[1]=="x" && colnames(sensfun)[2] == "var")
   Sens <- sensfun[,-(1:2)] else Sens<-sensfun

npar <- ncol(Sens)
L2 <- sqrt(colSums(Sens*Sens))

iNa <- 0
# Check for non-identifiable parameters
if (any(L2==0)){
  iNa<-which(L2==0)
  warning (paste("Sensitivity of parameter",colnames(Sens)[iNa],"is 0! "))
}
if (npar > 14)
    warning ("will reduce collinearity estimates: too many combinations")

normSens <- t(t(Sens) / L2)


Collin <- NULL

# internal function to generate collinearities
# for a given set of parameter combinations
collFun <- function(cc)
{
  for (i in 1:nrow(cc))
  {
   ii    <- cc[i,]
   S     <- normSens[,ii]
   Nident <- (iNa != 0 & iNa %in% ii)
   if (Nident) id <- Inf
   else id    <- 1/sqrt(min(eigen(t(S)%*%S)$value))
   psub     <- rep(0,npar)
   psub[ii] <- 1
   n        <- ncol(cc)
   Collin   <<- rbind(Collin,c(psub,n,id))
  }

}
if (is.null(parset)){

  combin <- function(n, v)  # combinations of n elements from a vector p (length (p) > = n)
  {
   if (n == 1)
     matrix(data=v, ncol=1)
   else if (n >= length(v))
     matrix(data=v, nrow=1)
   else rbind(cbind(v[1],combin(n-1,v[-1])),combin(n,v[-1]))
  }

 pset   <- 1:npar

 for (n in 2:npar)
 {
   numcomb <- choose(npar,n)
   if (numcomb < 5000)
   {
   cc  <-combin(n,pset)
   collFun(cc)
   }
 }
} else
{
if (! is.vector(parset))
  stop("'parset' should be a vector")
parset<-matrix(nr=1,parset)
collFun(parset)
}
Collin <- as.data.frame(Collin)

class(Collin) <- c("collin","data.frame")
names(Collin)<-c(colnames(Sens),"N","collinearity")
return(Collin)
}

plot.collin <- function(x,...)
{
nc <- ncol(x)
plot(x[,nc-1],x[,nc],main="Collinearity",
     xlab="Nr parameters",ylab="-",...)
}

print.collin <- function (x,...)
print(format(as.data.frame(unclass(x)),digits=2,scientific=FALSE,...))
