#' TDI.Escaramis
#'
#' @param data.long name of the dataset containing at least 3 columns (measurement, subject effect, rater effect).
#' @param y measurement column name.
#' @param id subject effect column name. Must be a factor.
#' @param met rater effect column name. Must be a factor.
#' @param rep replication effect column name. Must be a factor. When there are no replicas rep = NA
#'            Default is NA.
#' @param p A value or vector of the percentiles for the TDI.
#' @param choose.model logical asking if the minimum AIC model should be chosen no matter the value of `int`.
#' @param int if FALSE, a model with no interaction is fitted; if TRUE, a model with interaction is fitted.
#'            Default is FALSE. Only makes sense for data with repetition
#' @param tol tolerance
#'            Default is 10^(-8)
#' @param alpha confidence level for the TDI's upper bounds.
#'              Default is 0.05.
#'
#' @importFrom nlme lme pdBlocked pdIdent getVarCov
#' @importFrom stats AIC qchisq pnorm qnorm qt
#'
#' @noRd


TDI.Escaramis <- function(data.long, y, id, met, rep = NA, p,
                          choose.model = T, tol = 10^(-8), int = F, alpha = 0.05){

  if (!is.na(NA)) {
    m =  length(levels(data.long[, rep]))
  } else{
    rep <- NA
    m = 1
  }

  # Number of subjects
  ns <- length(unique(data.long$id))

  if(choose.model){
    if(!is.na(rep)){
      model.no.int <- lme(y~met, data = data.long, random = ~1|id, method = "REML")
      model.int <- lme(y~met, data = data.long,
                       random = list(id = pdBlocked(list(~1, pdIdent(form = ~-1+met)))),
                       method = "REML")
      if(AIC(model.no.int) <= AIC(model.int)){
        model <- model.no.int
        md <- abs(summary(model)$coefficients$fixed[2])
        sd <- sqrt(2)*model$sigma # Standard deviation of the difference
        df <- 2*ns*m - (ns + m - 1) # Degrees of freedom
      } else{
        model <- model.int
        md <- abs(summary(model)$coefficients$fixed[2])
        var.int <- getVarCov(model)[2,2]
        sd <- sqrt(2*(var.int + model$sigma^2)) # Standard deviation of the difference
        df <- 2*ns*(m - 1) # Degrees of freedom
      }
    } else{
      model <- lme(y~met, data = data.long, random = ~1|id, method = "REML")
      md <- abs(summary(model)$coefficients$fixed[2])
      sd <- sqrt(2)*model$sigma # Standard deviation of the difference
      df <- 2*ns*m - (ns + m - 1) # Degrees of freedom
    }
  } else if(int == F){
    model <- lme(y~met, data = data.long, random = ~1|id, method = "REML")
    md <- abs(summary(model)$coefficients$fixed[2])
    sd <- sqrt(2)*model$sigma # Standard deviation of the difference
    df <- 2*ns*m - (ns + m - 1) # Degrees of freedom
  } else{
    model <- lme(y~met, data = data.long, random = list(id = pdBlocked(list(~1,pdIdent(form = ~-1+met)))), method = "REML")
    md <- abs(summary(model)$coefficients$fixed[2])
    var.int <- getVarCov(model)[2,2]
    sd <- sqrt(2*(var.int + model$sigma^2)) # Standard deviation of the difference
    df <- 2*ns*(m - 1) # Degrees of freedom
  }

  # tdi
  nc <- md/sd
  tdi <- sd*sqrt(qchisq(p,1,ncp = (nc^2)))

  # Computing p1 using binary search algorithm
  low <- p
  high <- 1
  while (low <= high){
    mid <- (high + low)/2
    p.est <- pnorm(qnorm(mid)) - pnorm(-qnorm(mid)-2*md/sd) - p
    if ((p.est <= tol)&(p.est >= -tol)){
      break
    }
    if (p.est > tol){
      high <- mid - tol
    }
    if (p.est < (-tol)){
      low <- mid + tol
    }
  }

  # ub
  n <- 2*m*ns
  k <- suppressWarnings(1/sqrt(n)*qt(1-alpha, df, ncp = (qnorm(mid)*sqrt(n))))
  ub <- md + k*sd
  result <- c(tdi, ub)
  out <- list(result = result, fitted.model = model)
  return(out)
}
