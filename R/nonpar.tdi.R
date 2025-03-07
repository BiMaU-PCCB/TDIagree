#'
#' @importFrom stats quantile qnorm sd
#' @importFrom boot boot
#'
#' @noRd

# auxiliar function for NonParTDI
d.boot <- function(dataset, i, prob, type){
  return(quantile(dataset[i], prob, type=type))
}


# d: vector of differences
# p: percentile for the TDI
# alpha: confidence level for the TDI's upper bound
# type: an integer between 1 and 9 selecting one of
#         the nine quantile algorithms (to be passed to
#         'quantile' function)
# R: The number of bootstrap replicates.


NonParTDI <- function(d, p, alpha=0.05, type=8, R=1000){

  tdi <- quantile(abs(d), p, type=type)
  # upper bound TDI
  b <- boot(abs(d), d.boot, R = R, stype = "i", prob=p, type=type)
  tdi.boot <- b$t

  ub.tdi.p <- quantile(tdi.boot, 1-alpha, type=type)

  if (any(tdi.boot==0)){
    warning(paste(sum(tdi.boot==0), "zeros found in the bootstrap vector, computation of UB.NLOG with original scale"),
            immediate.=T)
    ub.tdi.norm.log <- tdi+qnorm(1-alpha)*sd(tdi.boot)
  }
  else{
    ln.tdi.boot <- log(tdi.boot)
    ub.tdi.norm.log <- exp(log(tdi)+qnorm(1-alpha)*sd(ln.tdi.boot))
  }

  res <- c(tdi, ub.tdi.p, ub.tdi.norm.log)
  names(res) <- c(paste0("TDI(", p*100, "%)"),
                  paste0("UB.TDI.P(", (1-alpha)*100, "%)"),
                  paste0("UB.TDI.NLOG(", (1-alpha)*100, "%)"))
  return(res)
}
