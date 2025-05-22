#' TDI.Escaramis
#'
#' @param data.long name of the dataset containing at least 3 columns (measurement, subject effect, rater effect).
#' @param y measurement column name.
#' @param id subject effect column name.
#' @param met rater effect column name.
#' @param rep replication effect column name.
#' @param p A value or vector of the probabilities for the TDI.
#' @param ub logical asking whether the UBs should be computed.
#'           Default is TRUE.
#' @param choose.model logical asking if the minimum AIC model should be chosen no matter the value of int.
#'                     Default is TRUE.
#' @param int if FALSE, a model with no interaction is fitted; if TRUE, a model with interaction is fitted.
#'            Default is FALSE. Only makes sense for data with repetition
#' @param add.es name of the columns in data that should be added to the model.
#'               Default is NULL.
#' @param tol tolerance.
#'            Default is 10^(-8)
#' @param alpha confidence level for inference on the TDI.
#'              Default is 0.05.
#'
#' @importFrom nlme lme pdBlocked pdIdent getVarCov
#' @importFrom stats AIC qchisq pnorm qnorm qt as.formula
#'
#' @noRd


TDI.Escaramis <- function(data.long, rep, p, ub = TRUE,
                          choose.model = TRUE, tol = 10^(-8), int = F, add.es = NULL, alpha = 0.05){

  if (!is.na(rep)) {
    m =  nlevels(data.long[, "rep"])
  } else{
    rep <- NA
    m = 1
  }

  ns <- length(unique(data.long$id))

  if (is.null(add.es)) {
    if (choose.model) {
      if (!is.na(rep)) {
        model.no.int <- lme(y~met, data = data.long, random = ~1|id, method = "REML")
        model.int <- lme(y~met, data = data.long,
                         random = list(id = pdBlocked(list(~1, pdIdent(form = ~-1+met)))),
                         method = "REML")
        if (AIC(model.no.int) <= AIC(model.int)) {
          model <- model.no.int
          md <- abs(summary(model)$coefficients$fixed[2])
          sd <- sqrt(2)*model$sigma
          df <- 2*ns*m - (ns + m - 1)
        } else{
          model <- model.int
          md <- abs(summary(model)$coefficients$fixed[2])
          var.int <- getVarCov(model)[2,2]
          sd <- sqrt(2*(var.int + model$sigma^2))
          df <- 2*ns*(m - 1)
        }
      } else{
        model <- lme(y~met, data = data.long, random = ~1|id, method = "REML")
        md <- abs(summary(model)$coefficients$fixed[2])
        sd <- sqrt(2)*model$sigma
        df <- 2*ns*m - (ns + m - 1)
      }
    } else if (int == F) {
      model <- lme(y~met, data = data.long, random = ~1|id, method = "REML")
      md <- abs(summary(model)$coefficients$fixed[2])
      sd <- sqrt(2)*model$sigma
      df <- 2*ns*m - (ns + m - 1)
    } else{
      model <- lme(y~met, data = data.long, random = list(id = pdBlocked(list(~1,pdIdent(form = ~-1+met)))), method = "REML")
      md <- abs(summary(model)$coefficients$fixed[2])
      var.int <- getVarCov(model)[2,2]
      sd <- sqrt(2*(var.int + model$sigma^2))
      df <- 2*ns*(m - 1)
    }
  } else{
    formula_str <- paste0("y ~ met + ", paste(add.es, collapse = " + "))
    formula_obj <- as.formula(formula_str)
    if (choose.model) {
      if (!is.na(rep)) {
        model.no.int <- lme(formula_obj, data = data.long, random = ~1|id, method = "REML")
        model.int <- lme(formula_obj, data = data.long,
                         random = list(id = pdBlocked(list(~1, pdIdent(form = ~-1+met)))),
                         method = "REML")
        if (AIC(model.no.int) <= AIC(model.int)) {
          model <- model.no.int
          md <- abs(summary(model)$coefficients$fixed[2])
          sd <- sqrt(2)*model$sigma
          df <- 2*ns*m - (ns + m - 1)
        } else{
          model <- model.int
          md <- abs(summary(model)$coefficients$fixed[2])
          var.int <- getVarCov(model)[2,2]
          sd <- sqrt(2*(var.int + model$sigma^2))
          df <- 2*ns*(m - 1)
        }
      } else{
        model <- lme(formula_obj, data = data.long, random = ~1|id, method = "REML")
        md <- abs(summary(model)$coefficients$fixed[2])
        sd <- sqrt(2)*model$sigma
        df <- 2*ns*m - (ns + m - 1)
      }
    } else if (int == F) {
      model <- lme(formula_obj, data = data.long, random = ~1|id, method = "REML")
      md <- abs(summary(model)$coefficients$fixed[2])
      sd <- sqrt(2)*model$sigma
      df <- 2*ns*m - (ns + m - 1)
    } else{
      model <- lme(formula_obj, data = data.long, random = list(id = pdBlocked(list(~1, pdIdent(form = ~-1+met)))), method = "REML")
      md <- abs(summary(model)$coefficients$fixed[2])
      var.int <- getVarCov(model)[2,2]
      sd <- sqrt(2*(var.int + model$sigma^2))
      df <- 2*ns*(m - 1)
    }
  }

  # tdi
  nc <- md/sd
  tdi.hat <- sd*sqrt(qchisq(p, 1, ncp = (nc^2)))

  ub.hat <- NA
  if(ub){
    # computing p1 using binary search algorithm
    low <- p
    high <- 1
    while (low <= high){
      mid <- (high + low)/2
      p.est <- pnorm(qnorm(mid)) - pnorm(-qnorm(mid)-2*md/sd) - p
      if ((p.est <= tol) & (p.est >= -tol)){
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
    ub.hat <- md + k*sd
  }
  result <- c(tdi.hat, ub.hat)
  out <- list(result = result, fitted.model = model)
  return(out)
}
