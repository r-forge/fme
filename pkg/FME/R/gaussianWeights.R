## This function adds a gaussian weighing function to a dataset in long format. By adding weights in a column called err. At the same time the dataset is
## extended to replicate each observation as many times as needed (see details).
##
## Parameters:
## obs         : data set in long format as is typically used in modCost
## x           : name of the independent variable in the dataset
## xmodel      : unique times at which model output is produced
## spread      : standard deviation used to calculate the weights from a normal distribution (stdev = standard deviation of that normal distribution)
## ordering    : list of variables in obs used for grouping during scale factor calculation (cf. tapply)
## weight      : scaling factor of modCost function (sd, mean, or none)

gaussianWeights <- function(obs,                            # data set in long format as is typically used in modCost (e.g. Meso1 in long format)
                            x           = "x",              # name of the independent variable (e.g. "day" in Meso1)
                            y           = "y",
                            xmodel      = NULL,             # unique times at which model output is produced
                            spread      = 1/24,             # standard deviation used to calculate the weights from a normal distribution (stdev = standard deviation of that normal distribution)
                            weight      = "none",           # scaling factor of modCost function (sd, mean, or none)
                            aggregation = "name",
                            ordering    = c(aggregation,x)  # list of variables in obs used for grouping during scale factor calculation (cf. tapply)
                           )
{

  # Exception/error handling
  if(is.null(xmodel))
    stop("Need times of model run to determine relevant time points")
  if(!(x %in% names(obs)))
    stop("Independent variable in observations not known")
  if(!(weight %in% c("mean","sd","none")))
    stop("Unknown weighing")

  # Selection of indices of relevant model times
  obsTimesSelect <- sapply(unique(obs[,x]),
                           FUN = function(t) which(abs(xmodel - t) == 0)
                           )
  obsRangeSelect <- sapply(unique(obs[,x]),
                           FUN = function(t) which(abs(xmodel - t) <= 3*spread)
                           )
  # Collection of actual model times
  obsTimes   <- lapply(obsTimesSelect, function(t) xmodel[t])
  obsRanges  <- lapply(obsRangeSelect, function(t) xmodel[t])
  obsWeights <- obsRanges # just to copy data structure

  # Calculation of inverse weights
  for(i in 1:length(obsTimes)) {
     obsWeights[[i]] <- dnorm(obsWeights[[i]],mean=obsTimes[[i]],sd=spread)
     obsWeights[[i]] <- obsWeights[[i]]/sum(obsWeights[[i]]) # rescaling to sum up to 1
     obsTimes[[i]] <- rep(obsTimes[[i]],times=length(obsWeights[[i]])) # replication of original observation times
  }

  # Construction of error data frame
  errors <- data.frame(time = unlist(obsTimes),
                       err = 1/unlist(obsWeights),
                       xmodel=unlist(obsRanges)
                      )
  # Merging with observations to obtain an extended observational dataset
  fullObs <- merge(obs,errors,by.x=x,by.y="time",all.y=T)
  fullObs[,x] <- fullObs[,"xmodel"]   # observational times are replaced by model times
  fullObs[,"xmodel"] <- NULL          # the latter column is redundant, hence removed

  if(weight != "none") {
    aggregationFunction <- tapply(fullObs[,y],INDEX=fullObs[,rev(aggregation)],weight,na.rm=T)
    factorLevels <- dimnames(aggregationFunction)
    names(factorLevels) <- rev(aggregation)
    aggregationFactor <- data.frame(expand.grid(factorLevels),
                                    scaling  = abs(array(aggregationFunction))
                                   )
    fullObs <- merge(fullObs,aggregationFactor)
    fullObs$err <- fullObs$err * fullObs$scaling
    fullObs$scaling <- NULL
  }

  fullObs <- eval(parse(text=paste("fullObs[order(",
                                   paste("fullObs$",ordering,sep="",collapse=","),
                                   "),]",
                                   collapse=""
                                  )
                       )
                 )

  firstColumns <- c(aggregation[1],x,y,"err")
  lastColumns  <- setdiff(names(fullObs),firstColumns)
  fullObs <- data.frame(fullObs[,firstColumns],fullObs[,lastColumns])
  names(fullObs) <- c(firstColumns,lastColumns)

  return(fullObs)
}
