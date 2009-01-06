
modFit <- function(f,p,...,lower=-Inf,upper=Inf,
                   method=c("Marq","Port","Newton","Nelder-Mead", "BFGS", "CG",
                   "L-BFGS-B", "SANN","Pseudo"),
                   control=list(),hessian=TRUE)
{ # check if valid input...
  method <- match.arg(method)

  pnames <- names(p)

  Lower <- rep(lower,len=length(p))
  Upper <- rep(upper,len=length(p))

# are boundaries allowed in the method?
  bounds <- method %in% c("L-BFGS-B", "Port", "Pseudo")

  FF      <- NULL
  useCost <- method != "Marq"   # marquardt uses residuals, others only model cost

  Func <- function(p,...)
   {
     FF<<- f(p,...)
     cM      <- class(FF) == "modCost"

     if (cM && useCost ) return(FF$model)
     if (cM)             return(FF$residual$res)
     if (useCost) return (sum(FF^2)) else return(FF)

   }
  Pars <- p

# Adapt function call if necessary
  if (bounds)                            # no need to change parameters
   Fun <- function(p,...) Func(p,...)
  else
    {
    # 1. Adapt initial parameter values
    lower <- -Inf
    upper <- Inf

      # lower and upper bounds...
    lu   <- which(is.finite(Lower) & is.finite(Upper))
    Pars[lu] <- tan(pi*((p[lu]-Lower[lu])/(Upper[lu]-Lower[lu]) - 0.5))
      # just lower bounds...
    l <- which(is.finite(Lower) & !is.finite(Upper))
    Pars[l] <-log(p[l]-Lower[l])
      # just upper bounds...
    u <- which(!is.finite(Lower) & is.finite(Upper))
    Pars[u] <- log(-p[u]+Upper[u])

    # 2. parameter transformation in function call
     Fun <- function(p,...)
     {
      PP     <- p
      PP[lu] <- Lower[lu]+(Upper[lu]-Lower[lu])*(atan(p[lu])/pi + 0.5)
      PP[l]  <- Lower[l]+exp(p[l])
      PP[u]  <- Upper[u]-exp(p[u])
      names(PP) <-pnames                #nlminb 'forgets' parameter names...
      Func(PP,...)
     }
    }

# optimise...
   if (method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))
   {
      res <- optim(Pars,Fun,method=method,lower=lower,upper=upper,
                   control=control,hessian=hessian,...)
      names(res)[2] <- "ssr"  #is called "value" here
      names(res)[3] <- "evaluations"  #is called "counts" here

   }
   else if (method == "Port")
   {
      res <- nlminb(start=Pars,objective=Fun,lower=lower,upper=upper,
                    control=control,...)
      names(res)[2] <- "ssr"  #is called "objective" here
      names(res)[6] <- "counts"
      # hessian not estimated...
      if(hessian) res$hessian <- hessian(Fun,res$par,centered=TRUE,...)
   }
   else if (method == "Marq")
   {
      Contr <- nls.lm.control()
      nmsC <- names(Contr)
      Contr[(namc <- names(control))]<-control
      if (length(noNms <- namc[!namc %in% nmsC]) > 0)
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

      res <- as.list(nls.lm(par=Pars,fn=Fun,control=Contr,...)[])

      names(res)[7] <- "iterations"  #is called "niter" here
      names(res)[9] <- "ssr"  #is called "deviance" here
      names(res)[3] <- "residuals"
      res$hessian <- 2*res$hessian  # returns 0.5*hessian!
      Diag<-unlist(res[6])
      res$diag<-NULL
      res$diag<-Diag
   }
   else if (method == "Newton")
   {
      res <- nlm(p=Pars,f=Fun,hessian=hessian,...)
      names(res)[1] <- "ssr"  #is called "minimum" here
      names(res)[2] <- "par"  #is called "estimate" here
   }
   else if (method == "Pseudo")
   {
      res <- pseudoOptim(p=Pars,f=Fun,lower=lower,upper=upper,
                         control=control,...)
      names(res)[2] <- "ssr"  #is called "cost" here
      if(hessian) res$hessian <- hessian(Fun,res$par,centered=TRUE,...)
   }
  if (! bounds && length(c(lu,l,u))>0)
  {
    respar <- res$par
    res$par[lu]<-Lower[lu]+(Upper[lu]-Lower[lu])*(atan(respar[lu])/pi + 0.5)
    res$par[l] <-Lower[l]+exp(respar[l])
    res$par[u] <-Upper[u]-exp(respar[u])
    # hessian re-estimated using backtransformed parameter values...
    useCost <- TRUE
    if(hessian) res$hessian <- hessian(Func,res$par,centered=TRUE,...)

  }
  names(res$par)<-names(p)
  class(res)    <- "modFit"
  if (!method == "Marq")
   if (class(FF) == "modCost") res$residuals <- FF$residuals$res else res$residuals  <- FF
  res
}

summary.modFit <- function (object, ...)  #inspired by summary.nls.lm
{
    param  <- object$par
    pnames <- names(param)
    p      <- length(param)
    covar  <- solve(0.5*object$hessian)   # unscaled covariance

    rdf <- length(object$residuals) - p
    resvar <- object$ssr / rdf
    se     <- sqrt(diag(covar) * resvar)
    names(se) <- pnames
    tval <- param/se

    param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(param) <- list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans <- list(residuals = object$residuals, variance=resvar, sigma = sqrt(resvar),
                df = c(p, rdf), cov.unscaled = covar,
                cov.scaled=covar * resvar,
                info = object$info, niter = object$iterations,
                stopmess = object$message,
                par = param)
    class(ans) <- "summary.modFit"
    ans
}
print.summary.modFit <-
  function (x, digits = max(3, getOption("digits") - 3), ...)

{
  df <- x$df
  rdf <- df[2]
  cat("\nParameters:\n")
  printCoefmat(x$par, digits = digits, ...)
  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")

  cat("\nParameter correlation:\n")
  printCoefmat(cov2cor(x$cov.unscaled), digits = digits, ...)

  invisible(x)
}
