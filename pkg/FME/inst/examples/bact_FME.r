
##################################################################
# A simple model of bacteria, growing on a substrate
# in a batch culture
# Model from Soetaert and Herman, 2009
##################################################################
mf <- par(mfrow=c(2,2))
require(FME)

##===============================================================##
##===============================================================##
##                           the Model                           ##
##===============================================================##
##===============================================================##

model <- function(t,state,pars)
{
with (as.list(c(state,pars)), {

dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact

return(list(c(dBact,dSub)))
                              })
}
# parameters
pars <- list(gmax =0.5,eff = 0.5,
              ks =0.5, rB =0.01, dB =0.01)

#initial conditions
Bini=0.1
Sini=100

# output times
tout    <- seq(0,50,by=0.5)
state   <- c(Bact=Bini,Sub = Sini)

# run model dynamically
out     <- as.data.frame(ode(state,tout,model,pars))

# plot output
plot(out$time,out$Bact,ylim=range(c(out$Bact,out$Sub)),
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)
lines(out$time,out$Sub,lty=2,lwd=2)
lines(out$time,out$Sub+out$Bact)

legend("topright",c("Bacteria","Glucose","TOC"),
       lty=c(1,2,1),lwd=c(2,2,1))

##===============================================================##
##===============================================================##
##       Global sensitivity analysis : Sensitivity ranges        ##
##===============================================================##
##===============================================================##

# the sensitivity parameters
parRanges <- data.frame(min=c(0.4,0.4,0.),max=c(0.6,0.6,0.02))
rownames(parRanges)<- c("gmax","eff","rB")
parRanges

# 1. Sensitivity parameter range: rB
#  equally-spaced parameters  "dist=grid"
tout    <- 0:50

# run model 100 times
Sens <- sensRange(func=model,y=state,times=tout,parms=pars,dist="grid",
                   sensvar="Bact",parRange=parRanges[3,],num=100)$summ

# plot output
yrange<-range(cbind(Sens$Min,Sens$Max))
plot(tout,ylim=yrange,xlab="time, hour",ylab="molC/m3",type="n",
     main="Sensitivity rB")

polygon(c(tout,rev(tout)),c(Sens$Min,rev(Sens$Max)),
        col=grey(0.9),border=NA)
polygon(c(tout,rev(tout)),c(Sens$Mean-Sens$Sd,
        rev(Sens$Mean+Sens$Sd)),col=grey(0.8),border=NA)
lines(tout,Sens$Mean,lwd=2)
legend("topleft",fill=c(grey(0.9),grey(0.8)),
       legend=c("Min-Max","Mean+-sd"),bty="n")
legend("left",lty=1,lwd=2,legend="Mean",bty="n")

# 2. sensitivity to all; latin hypercube sampling
Sens2 <- sensRange(func=model,y=state,times=tout,parms=pars,dist="latin",
                   sensvar="Bact",parRange=parRanges,num=100)$summ
yrange<-range(cbind(Sens2$Min,Sens2$Max))
plot(tout,ylim=yrange,xlab="time, hour",ylab="molC/m3",type="n",
     main="Sensitivity gmax,eff,rB")

polygon(c(tout,rev(tout)),c(Sens2$Min,rev(Sens2$Max)),
        col=grey(0.9),border=NA)
polygon(c(tout,rev(tout)),c(Sens2$Mean-Sens2$Sd,
          rev(Sens2$Mean+Sens2$Sd)),col=grey(0.8),border=NA)
lines(tout,Sens2$Mean,lwd=2)

par(mfrow=mf )

##===============================================================##
##===============================================================##
##       Local sensitivity analysis : sensitivity functions      ##
##===============================================================##
##===============================================================##

Sns<- sensFun(func=model,y=state,times=tout,parms=pars,
                   sensvar="Bact",varscale=1)

SnsBact <- Sns$fun

matplot(tout,SnsBact[,-(1:2)],type="l",ylab="Sfun",
        main="sensitivity functions")

Sns$model

SF<- sensFun(func=model,y=state,times=tout,parms=pars,
             sensvar=c("Bact","Sub"),varscale=1)
SF$model
SF$var
# Bivariate sensitivity

panel.cor <- function(x, y)
             text(x=mean(range(x)),y=mean(range(y)),
             labels=format(cor(x,y),digits=2))

pairs(SnsBact[,-(1:2)],upper.panel=panel.cor)
mtext(outer=TRUE,side=3,line=-2,
      "Sensitivity functions",cex=1.5)

# pairwise correlation
cor(SnsBact[,-(1:2)])

###############################
# multivariate sensitivity analysis
###############################
Coll <- collin(SnsBact[,-(1:2)])

# The larger the collinearity, the less identifiable the data set
format(data.frame(Coll),digits=2,scientific=FALSE)

nc <- ncol(Coll)
plot(Coll[,nc-1],Coll[,nc],log="y",main="Collinearity",
     xlab="Nr parameters",ylab="coll")

# 20 = magical number above which there are identifiability problems
abline(h=20,col="red")
Coll [Coll[,"collinearity"]<20&Coll[,"N"]==4,]


##===============================================================##
##===============================================================##
##         Fitting the model to the data - using nlminb          ##
##===============================================================##
##===============================================================##

# the "data"
Data <- matrix (nc=2,byrow=2,data=
c(  2,  0.14,    4,  0.21,    6,  0.31,    8,  0.40,
   10,  0.69,   12,  0.97,   14,  1.42,   16,  2.0,
   18,  3.0,    20,  4.5,    22,  6.5,    24,  9.5,
   26, 13.5,    28, 20.5,    30,  29 , 35, 65, 40, 61)
)
colnames(Data) <- c("time","Bact")
head(Data)

# 1. Objective function to minimise; parameters gmax and eff are fitted

Objective <- function (x)
{
 pars[c("gmax","eff")]<- x

 # output times
 tout    <- seq(0,50,by=0.5)
 state   <- c(Bact=Bini,Sub = Sini)

 # run model dynamically
 out     <- as.data.frame(ode(state,tout,model,pars))

 # Model cost
 Cost  <- modCost(obs=Data,model=out)

 return(Cost$model)
}

# 2. nlminb finds the minimum; parameters constrained to be > 0
print(system.time(Fit<-nlminb(start=c(0.5,0.5),
                  obj=Objective,lower=c(0.,0.))))

Fit

# Run best-fit model....
 pars[c("gmax","eff")]<- Fit$par

 # output times
 tout    <- seq(0,50,by=0.5)
 state   <- c(Bact=Bini,Sub = Sini)

 # run model dynamically
 out     <- as.data.frame(ode(state,tout,model,pars))

 # Model cost
 Cost  <- modCost(obs=Data,model=out)
Cost

# Plot residuals
plot(Cost$residual$x, Cost$residual$res,xlab="time",ylab="",main="residual")

# plot output
plot(out$time,out$Bact,ylim=range(out$Bact),
     xlab="time, hour",ylab="molC/m3",type="l",lwd=2)
points(Data,cex=2,pch=18)
