#' TDI.Choudhary.NP
#'
#' @param data name of the data set containing at least 3 columns (measurement, subject effect, rater effect).
#' @param cont measurement column name.
#' @param id subject effect column name.
#' @param met rater effect column name.
#' @param rep replication effect column name.
#'            Default is NA.
#' @param p A value or vector of the percentiles for the TDI.
#' @param alpha confidence level for the TDI's upper bounds.
#'              Default is 0.05.
#'
#' @importFrom stats confint
#' @importFrom multcomp glht parm
#'
#' @noRd


ChoudNonPar <- function(data, cont, id, met, rep = NA, p = 0.95, alpha = 0.05){

  d1 <- data[data[, met] == levels(data[, met])[1],]
  d2 <- data[data[, met] == levels(data[, met])[2],]

  d.12 <- drop(data.frame(cbind(subid=d1[, id], met1=d1[, cont], met2=d2[, cont])))

  if(is.na(rep)){
    d.paired <- with(d.12, data.frame(subid=subid, x=met1, y=met2))
    nrep <- 1
  }
  else{
    nrep <- length(levels(data[, rep]))

    d <- split(d.12, d.12$subid)
    d.paired <- data.frame(do.call("rbind", lapply(d, make.pairs)))
  }
  nsub <- length(unique(data[, id]))
  wt <- 1/(nsub*(nrep^2))

  # tdi
  emp.dist <- cbind(d.paired, jt.pmf=rep(wt, nsub*(nrep^2)))
  tdi <- tdi.fun(emp.dist[, 2:4], p) #

  # ub
  cdf.est <- cdf.fun(emp.dist[,2:4], tdi) #
  infl.cdf <- cdf.infl.fun(emp.dist[,2:4], p) #
  cdf.var <- asymp.var(emp.dist, infl.cdf, nsub = nsub, nrep = nrep, wt = wt)/nsub #

  cdf.glht <- confint(glht(model = parm(coef = cdf.est, vcov = as.matrix(cdf.var)),
                           linfct = diag(1), alternative = "less"), level = (1-alpha))
  cval.cdf <- attr(cdf.glht$confint, "calpha")

  ub <- as.numeric(tdi.fun(emp.dist[, 2:4], (p+cval.cdf*sqrt(cdf.var))))
  out <- c(tdi, ub)
  return(out)

}
