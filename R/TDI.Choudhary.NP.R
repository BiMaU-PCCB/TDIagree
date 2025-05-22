#' TDI.Choudhary.NP
#'
#' @param data name of the data set containing at least 3 columns (measurement, subject effect, rater effect).
#' @param y measurement column name.
#' @param id subject effect column name.
#' @param met rater effect column name.
#' @param rep replication effect column name.
#' @param p a value or vector of the probabilities for the TDI.
#' @param ub logical asking whether the UBs should be computed.
#'           Default is TRUE.
#' @param alpha confidence level for inference on the TDI.
#'              Default is 0.05.
#'
#' @importFrom stats aggregate confint
#' @importFrom multcomp glht parm
#'
#' @noRd


# auxiliar functions for TDI.Choudhary.NP

make.pairs <- function(d) {
  data <- as.data.frame(d)
  nrep <- dim(data)[1]
  stopifnot(nrep > 1)
  d.paired <- cbind(id = rep(data[, 1], nrep), y_met1 = rep(data[, 2], each = nrep),
                    y_met2 = rep(data[, 3], nrep))
  d.paired
}

square.of.sum <- function(x, nrep) {
  x <- matrix(x, nrow = nrep, ncol = nrep, byrow = T)
  rsum <- sum(rowSums(x)^2)
  csum <- sum(colSums(x)^2)
  asum <- sum(x)^2
  c(rsum, csum, asum)
}



TDI.Choudhary.NP <- function(data, rep, p, ub = TRUE, alpha = 0.05){

  d1 <- data[data[, "met"] == levels(data[, "met"])[1], ]
  d2 <- data[data[, "met"] == levels(data[, "met"])[2], ]

  d.12 <- drop(data.frame(cbind(id = d1[, "id"], y_met1 = d1[, "y"], y_met2 = d2[, "y"])))
  if (is.na(rep)) {
    d.paired <- with(d.12, data.frame(id = id, y_met1 = y_met1, y_met2 = y_met2))
    nrep <- 1
  } else{
    nrep <- nlevels(data[, "rep"])

    d <- split(d.12, d.12$id)
    d.paired <- data.frame(do.call("rbind", lapply(d, make.pairs)))
  }
  nsub <- length(unique(data[, "id"]))
  wt <- 1/(nsub*(nrep^2))

  # tdi
  emp.dist <- cbind(d.paired, jt.pmf = rep(wt, nsub*(nrep^2)))

  ad <- abs(emp.dist$y_met1 - emp.dist$y_met2)
  ad.pmf <- aggregate(emp.dist$jt.pmf, by = list(ad = ad), FUN = sum)
  ad.cdf <- cumsum(ad.pmf[, 2])
  if (isTRUE(all.equal(max(ad.cdf), 1))) {ad.cdf[length(ad.cdf)] <- 1 }

  tdi.hat <- min(ad.pmf[ad.cdf >= p, 1])

  # ub
  ub.hat <- NA
  if (ub) {
    cdf.est <- sum(emp.dist$jt.pmf[(ad <= tdi.hat)])
    infl.cdf <- cdf.est - 1*(ad <= tdi.hat)

    ssq <-  sum(infl.cdf^2)
    ss <- rowSums(sapply(split(infl.cdf, emp.dist$id), FUN = "square.of.sum", nrep = nrep))
    # [1] = row, [2] = col, [3] = all
    if (nrep != 1){
      variance <- wt*ssq
      cov1 <- (ss[1] - ssq)/(nsub*(nrep^2)*(nrep-1))
      cov2 <- (ss[2] - ssq)/(nsub*(nrep^2)*(nrep-1))
      cov3 <- (ss[3] - ss[1] - ss[2] + ssq)/(nsub*(nrep^2)*((nrep-1)^2))
      cdf.var <- (variance + (nrep-1)*cov1 + (nrep-1)*cov2 + ((nrep-1)^2)*cov3)*wt
    } else{
      variance <- wt*ssq
      cdf.var <- variance/nsub
    }

    cdf.glht <- confint(glht(model = parm(coef = cdf.est, vcov = as.matrix(cdf.var)),
                             linfct = diag(1), alternative = "less"), level = (1-alpha))
    cval.cdf <- attr(cdf.glht$confint, "calpha")

    c <- p+cval.cdf*sqrt(cdf.var)
    if (c > 1) {c <- 1}
    ub.hat <- min(ad.pmf[ad.cdf >= c, 1])
  }
  out <- c(tdi.hat, ub.hat)
  return(out)
}
