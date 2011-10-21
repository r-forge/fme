\name{gaussianWeights}
\alias{gaussianWeights}

\title{
  A kernel average smoother function to weigh residuals according to a Gaussian
  density function
}

\description{
  A calibration dataset in database format (cf. modCost for the database
  format) is extended in order to fit model output surrounding a specific
  observation in the calibration dataset to this observation by weighing
  the residuals according to a Gaussian density function that has the
  relevant observation as mean and a specified standard deviation.
  
  This weighted regression approach may be interesting when discontinuities
  exist in the observational data. Under these circumstances small changes 
  in the timing (or more general the position along the axis of the independent
  variable) of the model output may have a disproportional impact on the
  overall goodness-of-fit (e.g. timing of nutrient depletion). Another
  application is when the exact value of the independent variable in the
  observational dataset is uncertain (e.g. slices of sediment profiles,
  or the timing of a sampling).

}
\usage{
gaussianWeights (obs, x = "x", y = "y", xmodel = NULL, spread = 1/24,
                weight = "none",aggregation = "name",ordering = c(aggregation,x))
}

\arguments{

  \item{obs }{dataset in long (database) format as is typically used in
     modCost
  }
  \item{x }{name of the independent variable (typically x, cf. modCost) in
     \code{obs}
  }
  \item{y }{name of the dependent variable in \code{obs}
  }
  \item{xmodel }{an ordered vector of unique times at which model output
     is produced
  }
  \item{spread }{standard deviation used to calculate the weights from a
     normal density function. This value also determines the number of
     point from the model output that are compared to a specific observa-
     tion in \code{obs} (2 * 3 * spread + 1; containing 99.7\% of the
     Gaussian distribution, centered around the observation of interest).
  }
  \item{weight }{scaling factor of the modCost function ("sd", "mean", or
    "none"). The Gaussian weights are multiplied by this factor.
  }
  \item{aggregation }{vector of column names from the dataset that are used
    to aggregate the scaling factor
  }
  \item{ordering }{list of variables in obs used for grouping during scale
    factor calculation (cf. tapply)
  }
}
\value{
  
  A modified version of \code{obs} is returned with the following extensions:

  1. Each observation is replicated n times were n represents the number
  of \code{modelx} values that are within 3 * the standard deviation of 
  the Gaussian density function with the \code{x} value of the relevant
  observation as mean and \code{spread} as standard deviation.

  2. These replicate observations get the same \code{x} values as their
  modeled counterparts.

  3. Weights are given in column, called "err"

  The returned data frame has the following columns:
  \itemize{
    \item{"name" or another name specified by the first element of 
    \code{aggregation}. Usually this column contains the names of the
    observed variables.
    }
    \item{"x" or another name specified by \code{x}
    }
    \item{"y" or another name specified by \code{y}
    }
    \item{"err" containing the calculated weights
    }
    \item{The rest of the columns of the data frame given by \code{obs} in
    that order
    }
  }
}

\details{

  !formules hier geven!

}

\examples{
## =======================================================================
## A Sediment example
## =======================================================================

## Suppose we have sediment oxygen profiles with measurements every
## centimeter in 3 sediment types
depth <- 0:7
observations <- data.frame(
                    profile = rep(c("mud","silt","sand"), each=8),
                    depth   = depth,
                    O2      = c(c(6,1,0.5,0.1,0.05,0,0,0),
                                c(6,5,3,2,1.5,1,0.5,0),
                                c(6,6,5,4,3,2,1,0)
                              )
                )

## The model generates profiles with a depth resolution of 1 millimeter
modeldepths <- seq(0,9,by=0.05)

gaussianWeights(obs=observations,x="depth",y="O2",xmodel = modeldepths,
                spread=0.1,weight="none",aggregation="profile")

}
\author{
  Tom Van Engeland <t.vanengeland@nioo.knaw.nl>
}

\keyword{utilities}

