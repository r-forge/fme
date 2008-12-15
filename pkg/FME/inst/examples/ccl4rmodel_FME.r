##################################################################
##################################################################
######            Applications of the ccl4 model            ######
##################################################################
##################################################################

par(mfrow=c(2,2))
require(FME)
##===============================================================##
##===============================================================##
##                           the Model                           ##
##===============================================================##
##===============================================================##

# the model is a FORTRAN-implemented example from deSolve

##====================================
## parameter values
##====================================

  Pm <- c(

   ### Physiological parameters
   BW= 0.182,   # Body weight (kg)
   QP= 4.0  ,   # Alveolar ventilation rate (hr^-1)
   QC= 4.0  ,   # Cardiac output (hr^-1)
   VFC= 0.08,   # Fraction fat tissue (kg/(kg/BW))
   VLC= 0.04,   # Fraction liver tissue (kg/(kg/BW))
   VMC= 0.74,   # Fraction of muscle tissue (kg/(kg/BW))
   QFC= 0.05,   # Fractional blood flow to fat ((hr^-1)/QC
   QLC= 0.15,   # Fractional blood flow to liver ((hr^-1)/QC)
   QMC= 0.32,   # Fractional blood flow to muscle ((hr^-1)/QC)

   ## Chemical specific parameters for chemical
   PLA= 16.17,  # Liver/air partition coefficient
   PFA= 281.48, # Fat/air partition coefficient
   PMA= 13.3,   # Muscle/air partition coefficient
   PTA= 16.17,  # Viscera/air partition coefficient
   PB= 5.487,   # Blood/air partition coefficient
   MW= 153.8,   # Molecular weight (g/mol)
   VMAX= 0.04321671, # Max. velocity of metabolism (mg/hr) -calibrated
   KM= 0.4027255,    # Michaelis-Menten constant (mg/l) -calibrated

   # Parameters for simulated experiment
   CONC= 1000,  # Inhaled concentration
   KL= 0.02,    # Loss rate from empty chamber /hr
   RATS= 1.0,   # Number of rats enclosed in chamber
   VCHC= 3.8    # Volume of closed chamber (l)
   )

##====================================
## state variables
##====================================
  y <- c(AI = 21,   # total mass , mg
         AAM = 0,
         AT = 0,
         AF = 0,
         AL = 0,
         CLT = 0, ### area under the conc.-time curve in the liver
         AM = 0   ### the amount metabolized (AM)
         )

##====================================
## the observations
##====================================

Obs <- ccl4data[ccl4data$initconc==1000,]
plot(ChamberConc ~ time,data=Obs,xlab="Time (hours)",
         xlim=range(c(0,ccl4data$time)),
         ylab="Chamber Concentration (ppm)",
         log="y",main = "ccl4model")

Pm["CONC"] <-1000

VCH <- Pm[["VCHC"]] - Pm[["RATS"]]*Pm[["BW"]]
AI0 <- VCH * Pm[["CONC"]]*Pm[["MW"]]/24450
y["AI"] <- AI0

##====================================
## run the model:
##====================================
times<-unique(Obs$time)
out <- as.data.frame(ccl4model(times,y,Pm))
lines(out$time,out$CP,lwd=2)

##===============================================================##
##===============================================================##
##       Local sensitivity analysis : sensitivity functions      ##
##===============================================================##
##===============================================================##

# 1. Sensitivity functions
# All parameters are sensitivity parameters, all variables selected
Sens <- sensFun(solver=ccl4model,times=times,y=y,parms=Pm,varscale=1)
Sens$model

##====================================
## univariate sensitivity
##====================================
format(Sens$model,digits=2)

# select the ones with highest sensitivity
pselect <- names(Pm)[which (Sens$model$L2>50)]

##====================================
## bivariate sensitivity
##====================================
panel.cor <- function(x, y)
             text(x=mean(range(x)),y=mean(range(y)),
             labels=format(cor(x,y),digits=2))

pairs(Sens$fun[,pselect],upper.panel=panel.cor,gap=0)
mtext(outer=TRUE,side=3,line=-1.5,
      "Sensitivity functions",cex=1.5)

##====================================
## multivariate sensitivity
##====================================

Coll<-collin(Sens$fun[,pselect])
head(Coll)
tail(Coll)

plot(Coll$N,Coll$collinearity,xlab="N",ylab="collinearity",
     main="identifiability",log="y")
abline(h=20,col="red")


# 'identifiable parameter combinations' with 7 parameters
Coll[which(Coll$N==7 & Coll$collinearity < 20),]

##===============================================================##
##===============================================================##
##       Global sensitivity analysis : Sensitivity ranges        ##
##===============================================================##
##===============================================================##

# 1. Define parameter range 0.8 - 1.2 their value (sensitive pars only)
parRange <- data.frame(min=Pm[pselect]*0.8,max=Pm[pselect]*1.2)
rownames(parRange) <- names(Pm[pselect])

parRange
# sensitivity range for sensitivity variable CP
Sr <- sensRange(solver=ccl4model,times=times,y=y,parms=Pm,sensvar="CP",
                parRange=parRange[1,],num=100)$summ

yrange<-range(cbind(Sr$Min,Sr$Max))
plot(times,xlim=range(times),ylim=yrange,xlab="time, hour",type="n",
     ylab="Chamber Concentration (ppm)",main="Sensitivity BW")

polygon(c(times,rev(times)),c(Sr$Min,rev(Sr$Max)),
        col=grey(0.9),border=NA)
polygon(c(times,rev(times)),c(Sr$Mean-Sr$Sd,
        rev(Sr$Mean+Sr$Sd)),col=grey(0.8),border=NA)
lines(times,Sr$Mean,lwd=2)
legend("topright",fill=c(grey(0.9),grey(0.8)),
       legend=c("Min-Max","Mean+-sd"),bty="n")
legend("right",lty=1,lwd=2,legend="Mean",bty="n")

# and so on for other parameters


##===============================================================##
##===============================================================##
##         Fitting the model to the data - using nlminb          ##
##===============================================================##
##===============================================================##

# fitted parameters
fitPAR <- c("VMAX","CONC","KM")

# model cost is to be minimised
Objective <- function(P)
{
 Pm[fitPAR]<-P
 VCH <- Pm[["VCHC"]] - Pm[["RATS"]]*Pm[["BW"]]
 AI0 <- VCH * Pm[["CONC"]]*Pm[["MW"]]/24450
 y["AI"] <- AI0

 out <- as.data.frame(ccl4model(times,y,Pm))
 Cost <- sum((out$CP-Obs$ChamberConc)^2)
 return(Cost)
}
# this one does not work
(Fit<-nlminb(start=c(VMAX=0.04,CONC=1000,KM=0.4),
             lower=c(0,0,0),obj=Objective))


##===============================================================##
##===============================================================##
## Fit the model to data using the Levenberg-Marquardt algorithm ##
##===============================================================##
##===============================================================##

# 1. Define the model residuals

Residuals <- function(P)
{
 Pm[fitPAR]<-P
 VCH <- Pm[["VCHC"]] - Pm[["RATS"]]*Pm[["BW"]]
 AI0 <- VCH * Pm[["CONC"]]*Pm[["MW"]]/24450
 y["AI"] <- AI0

 out <- as.data.frame(ccl4model(times,y,Pm))
 return(out$CP-Obs$ChamberConc)
}

(MrqFit<-nls.lm(par=c(VMAX=1,CONC=500,KM=1),fn=Residuals))
MrqFit$hessian  #!!!

# run model with the optimized value:
 Pm[fitPAR]<-MrqFit$par
 VCH <- Pm[["VCHC"]] - Pm[["RATS"]]*Pm[["BW"]]
 AI0 <- VCH * Pm[["CONC"]]*Pm[["MW"]]/24450
 y["AI"] <- AI0

 fitted <- as.data.frame(ccl4model(times,y,Pm))

plot(ChamberConc ~ time,data=Obs,xlab="Time (hours)",
         xlim=range(c(0,ccl4data$time)),
         ylab="Chamber Concentration (ppm)",
         log="y",main = "ccl4model-fitted")
lines(fitted$time,fitted$CP,lwd=2)

MrqFit[]
