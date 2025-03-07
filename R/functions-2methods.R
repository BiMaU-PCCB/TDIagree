# functions needed to compute tdi nonpar by Choudhary

#' @importFrom stats aggregate
#'

############
make.pairs <- function(idat)  # forms pairs of (x, y), idat has 3 cols (subid, x, y)
{
  idat <- as.data.frame(idat)
  irep <- dim(idat)[1]
  stopifnot(irep > 1)
  paired.dat <- cbind(subid=rep(idat[,1], irep), x=rep(idat[,2], each=irep),
                      y=rep(idat[,3], irep))
  paired.dat
}
################
ssquares.for.cov <- function(x12, x13, x23, nrep=3)
{
  x12 <- matrix(x12, nrow=nrep, ncol=nrep, byrow=T)
  x13 <- matrix(x13, nrow=nrep, ncol=nrep, byrow=T)
  x23 <- matrix(x23, nrow=nrep, ncol=nrep, byrow=T)

  sum12.13 <- sum(rowSums(x12)*rowSums(x13))
  sum12.23 <- sum(colSums(x12)*rowSums(x23))
  sum13.23 <- sum(colSums(x13)*colSums(x23))
  asum12.13 <- sum(x12)*sum(x13)
  asum12.23 <- sum(x12)*sum(x23)
  asum13.23 <- sum(x13)*sum(x23)

  c(sum12.13=sum12.13, sum12.23=sum12.23, sum13.23=sum13.23, asum12.13= asum12.13,
    asum12.23=asum12.23, asum13.23=asum13.23)
}
############
square.of.sum <- function(x, nrep) {
  x <- matrix(x, nrow=nrep, ncol=nrep, byrow=T)
  rsum <- sum(rowSums(x)^2)
  csum <- sum(colSums(x)^2)
  asum <- sum(x)^2
  c(rsum, csum, asum)
}
##############
tdi.fun <- function(dist, p) # dist has 3 cols - x, y, jt.pmf
  # inverse of cdf, agrees with quantile
{
  #if(p>1) {p <- 1}
  x1 <- dist[,1]
  x2 <- dist[,2]
  jt.pmf <- dist[,3]
  # dist of |d|
  ad <- abs(x1-x2)
  ad.pmf <- aggregate(jt.pmf, by=list(ad=ad), FUN=sum)  #col1 = ad, col2=pmf
  ad.cdf <- cumsum(ad.pmf[,2])
  if(isTRUE(all.equal(max(ad.cdf), 1))) {ad.cdf[length(ad.cdf)] <- 1 } # fix for max not exactly 1
  # tdi
  tdi <- min(ad.pmf[ad.cdf >= p, 1])
  tdi
}
##############
cdf.fun <- function (dist, u)  # dist has 3 columns: x1, x2, jt.pmf
{
  x1 <-  dist[,1]
  x2 <- dist[,2]
  jt.pmf <- dist[,3]
  ad <- abs(x1-x2)
  cdf <- sum(jt.pmf[(ad <= u)])
  cdf
}
##################
asymp.cov <- function(dist12, dist13, dist23, infl12, infl13, infl23, nsub, nrep) # each dist has 4 cols - subid, x, y, jt.pmf
{

  infl12 <- split(infl12, dist12$subid)
  infl13 <- split(infl13, dist13$subid)
  infl23 <- split(infl23, dist23$subid)

  ss <- rowSums(mapply(FUN="ssquares.for.cov", x12=infl12, x13=infl13, x23=infl23))

  # 12 vs 13
  cov1 <- ss["sum12.13"]/(nsub*(nrep^3))
  cov2 <- (ss["asum12.13"]-ss["sum12.13"])/(nsub*(nrep^3)*(nrep-1))
  cov12.13 <- ((nrep^3)*cov1 + (nrep^3)*(nrep-1)*cov2)/(nrep^4)


  # 12 vs 23
  cov1 <- ss["sum12.23"]/(nsub*(nrep^3))
  cov2 <- (ss["asum12.23"]-ss["sum12.23"])/(nsub*(nrep^3)*(nrep-1))
  cov12.23 <- ((nrep^3)*cov1 + (nrep^3)*(nrep-1)*cov2)/(nrep^4)

  # 13 vs 23
  cov1 <- ss["sum13.23"]/(nsub*(nrep^3))
  cov2 <- (ss["asum13.23"]-ss["sum13.23"])/(nsub*(nrep^3)*(nrep-1))
  cov13.23 <- ((nrep^3)*cov1 + (nrep^3)*(nrep-1)*cov2)/(nrep^4)

  result <- c(cov12.13, cov12.23, cov13.23)
  names(result) <- c("12.13", "12.23", "13.23")
  result
}
#################
# afegim nsub, nrep i wt com a arguments sinó s'han de declarar globalment
asymp.var <- function(dist, infl, nsub, nrep, wt) # dist has 4 cols - subid, x, y, jt.pmf
{
  ssq <-  sum(infl^2)
  ss <- rowSums(sapply(split(infl, dist$subid), FUN="square.of.sum", nrep=nrep))
  # [1] = row, [2] = col, [3] = all
  if (nrep!=1){
    ivar <- wt*ssq      # wt = 1/(nsub*nrep^2)
    icov1 <- (ss[1] - ssq)/(nsub*(nrep^2)*(nrep-1))
    icov2 <- (ss[2] - ssq)/(nsub*(nrep^2)*(nrep-1))
    icov3 <- (ss[3] - ss[1] - ss[2] + ssq)/(nsub*(nrep^2)*((nrep-1)^2))
    avar <- (ivar + (nrep-1)*icov1 + (nrep-1)*icov2 + ((nrep-1)^2)*icov3)/(nrep^2)
    avar
  }else{
    ivar <- wt*ssq
    avar <- ivar
    avar
  }
}
###############
cdf.infl.fun <- function(dist, p) # dist has 3 cols - x, y, jt.pmf
{
  x1 <-  dist[,1]
  x2 <- dist[,2]
  jt.pmf <- dist[,3]
  ad <- abs(x1-x2)
  tdi.hat <- tdi.fun(dist, p)
  infl <- cdf.fun(dist, tdi.hat) - 1*(ad <= tdi.hat)
  infl
}
################
make.pairs <- function(idat)  # forms pairs of (x, y), idat has 3 cols (subid, x, y)
{
  idat <- as.data.frame(idat)
  irep <- dim(idat)[1]
  stopifnot(irep > 1)
  paired.dat <- cbind(subid=rep(idat[,1], irep), x=rep(idat[,2], each=irep),
                      y=rep(idat[,3], irep))
  paired.dat
}
################
