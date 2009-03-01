
## -----------------------------------------------------------------------------
## Latin Hypercube distribution
## -----------------------------------------------------------------------------

Latinhyper <- function(parRange, num) {

  npar <- nrow(parRange)
  latin <- NULL
  for (i in 1:npar) {
    pr    <- unlist(parRange[i,])
    dpar  <- diff(pr)/num                             # delt of parameter interval
    ii    <- sort(runif(num),index.return=TRUE)$ix-1  # index to interval (0=1st interval)
    rr    <- runif(num)
    pval  <- pr[1]+ii*dpar+rr*dpar                    # random value within interval
    latin <- cbind(latin,pval)
  }
  colnames(latin) <- rownames(parRange)
  latin
}

