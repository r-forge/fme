## Function to convert a crosstab dataset to long format (database).
##
## Parameters:
## . data      : the crosstab dataset
## . x         : the independent variable to be replicated
## . include   : a character vector containing the names of the columns that should be included
## . exclude   : a character vector containing the names of the columns that should be excluded
## . replicate : a character vector containing the names of columns that have to be replicated for every column name (e.g. variable name)
## . error     : boolean indicating whether the final dataset in long format should contain an extra column for error values (standard deviation, standard error, ...)

cross2long <- function(data,
                       x,
                       include=NULL,
                       exclude=NULL,
                       replicate=NULL,
                       error=F
                      )
{
  # Exception/error handling
  data <- as.data.frame(data)
  nbUnique <- length(unique(c(include,exclude,replicate)))
  nbTotal  <- (length(include) + length(exclude) + length(replicate))
  if(nbUnique != nbTotal)
    stop("Include, exclude and replicate should be disjoint sets of column names")
  if(is.null(x))
    stop("An independent variable should be given")
  indep <- data[,x]
  data[x] <- NULL
  if(!is.null(replicate)) {
    toBeReplicated <- subset(data,select=replicate)
  }
  if(!is.null(exclude))
    data <- data[!(names(data) %in% exclude)]
  if(!is.null(include))
    data <- data[include]

  # gathering of data
  stackedData <- stack(data)
  nNames <- length(unique(stackedData$ind))
  nObs   <- nrow(stackedData)
  indep <- rep(indep,times=nNames)
  if(!is.null(replicate)) {
    replicated <- as.data.frame(lapply(as.list(toBeReplicated), FUN=function(x) rep(x,times=nNames)))
    # Construction of appropriate data structure
    if(error)
      out <- data.frame(name = factor(stackedData$ind,exclude=NULL),
                        x    = indep,
                        y    = stackedData$values,
                        err  = rep(1,nObs),
                        as.data.frame(lapply(replicated,function(x) factor(x,exclude=NULL)))
                       )
    else
      out <- data.frame(name = factor(stackedData$ind,exclude=NULL),
                        x    = indep,
                        y    = stackedData$values,
                        as.data.frame(lapply(replicated,function(x) factor(x,exclude=NULL)))
                       )
  }
  else {
    if(error)
      out <- data.frame(name = factor(stackedData$ind,exclude=NULL),
                        x    = indep,
                        y    = stackedData$values,
                        err  = rep(1,nObs)
                       )
    else
      out <- data.frame(name = factor(stackedData$ind,exclude=NULL),
                        x    = indep,
                        y    = stackedData$values
                       )

  }

  return(out)
}
