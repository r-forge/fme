\name{plotmod}
\alias{plotmod}
\alias{histmod}
\title{
  plots model output in separate plots
  }
\description{
  Produces x-y plots or histograms of the columns from a matrix in separate figures
}

\usage{
plotmod(mat, x=1, which = 2:ncol(mat), trace=FALSE, \dots)

histmod(mat, which = 2:ncol(mat), \dots)
}

\arguments{
  \item{mat }{Matrix whose columns have to be plotted, versus the \code{x}-column
  }
  \item{x }{the name or the index to the column of \code{mat}, which should
    be used as x-values. Default = first column.
  }
  \item{which }{the name or the index to the columns that should be plotted as
    y-values. Default = all columns except the first.
  }
  \item{trace }{if \code{TRUE}, adds smoothed line to the plot.
  }
  \item{...}{additional arguments passed to function \code{plot}.
  }
}

\details{
  The method will present the figures in n rows and n columns, except when
  \code{mfrow} is passed as extra argument

  Do not use with too many columns
}

\examples{

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plotting distributions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
n<-200
mat<-data.frame(cbind(iter=1:n, norm=rnorm(n), unif=runif(n),
  lnorm=rlnorm(n), logistic=rlogis(n)))
plotmod(mat)

histmod(mat,breaks=20,col="darkblue")
histmod(mat,which=c("norm","lnorm"),mfrow=c(1,2))

}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\keyword{utilities}

