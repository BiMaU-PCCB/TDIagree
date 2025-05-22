#' TDI.Choudhary.P
#'
#' @param data.long name of the dataset containing at least 3 columns (measurement, subject effect, rater effect).
#' @param y measurement column name.
#' @param id subject effect column name.
#' @param met rater effect column name.
#' @param rep replication effect column name.
#' @param p A value or vector of the probabilities for the TDI.
#' @param ub logical asking whether the UBs should be computed.
#'           Default is TRUE.
#' @param choose.model logical asking if the minimum AIC model should be chosen no matter the value of var.equal.
#'                     Default is TRUE.
#' @param var.equal logical asking if there is homoscedasticity between methods.
#'                  Default is TRUE.
#' @param alpha confidence level for inference on the TDI.
#'              Default is 0.05.
#'
#' @importFrom stats qchisq dnorm qnorm AIC
#' @importFrom nlme lme varIdent lmeControl
#'
#' @noRd


# auxiliar functions for TDI.Choudhary.P

TDI.Choudhary.P.EqVar <- function(p, mu1, mu2, sigma, invI, ub = ub, alpha = 0.05){
  # tdi
  mu.hat <- mu1 - mu2
  sigma.hat <- sqrt(2)*sigma
  tdi.hat <- sigma.hat*sqrt(qchisq(p, df = 1, ncp = (mu.hat/sigma.hat)^2))
  # ub
  ub.hat <- NA
  if (ub) {
    z.l <- -(tdi.hat + mu.hat)/sigma.hat
    z.u <- (tdi.hat - mu.hat)/sigma.hat
    s <- dnorm(z.l) + dnorm(z.u)
    partial.mu1 <- (1/(s*tdi.hat))*(dnorm(z.u) - dnorm(z.l))
    partial.mu2 <- (1/(s*tdi.hat))*(dnorm(z.l) - dnorm(z.u))
    partial.logsigma <- sigma.hat^2/(s*tdi.hat*sigma.hat)*(z.u*dnorm(z.u)-z.l*dnorm(z.l))
    G <- c(partial.mu1, partial.mu2, partial.logsigma)
    ub.hat <- exp(log(tdi.hat) + qnorm(1-alpha)*sqrt(t(G)%*%as.matrix(invI)%*%G))
  }
  res <- c(tdi.hat, ub.hat)
  return(res)
}

TDI.Choudhary.P.UneqVar <- function(p, mu1, mu2, sigma1, sigma2, invI, ub = ub, alpha = 0.05){
  # tdi
  mu.hat <- mu1 - mu2
  sigma.hat <- sqrt(sigma1^2 + sigma2^2)
  tdi.hat <- sigma.hat*sqrt(qchisq(p, df = 1, ncp = (mu.hat/sigma.hat)^2))
  # ub
  ub.hat <- NA
  if (ub) {
    z.l <- -(tdi.hat + mu.hat)/sigma.hat
    z.u <- (tdi.hat - mu.hat)/sigma.hat
    s <- dnorm(z.l) + dnorm(z.u)
    partial.mu1 <- (1/(s*tdi.hat))*(dnorm(z.u) - dnorm(z.l))
    partial.mu2 <- (1/(s*tdi.hat))*(dnorm(z.l) - dnorm(z.u))
    partial.logsigma1 <- sigma1^2/(s*tdi.hat*sigma.hat)*(z.u*dnorm(z.u)-z.l*dnorm(z.l))
    partial.logsigma2 <- sigma2^2/(s*tdi.hat*sigma.hat)*(z.u*dnorm(z.u)-z.l*dnorm(z.l))
    G <- c(partial.mu1, partial.mu2, partial.logsigma1, partial.logsigma2)
    ub.hat <- exp(log(tdi.hat) + qnorm(1-alpha)*sqrt(t(G)%*%as.matrix(invI)%*%G))
  }
  res <- c(tdi.hat, ub.hat)
  return(res)
}



TDI.Choudhary.P <- function(data.long, p, ub = TRUE,
                            choose.model = TRUE, var.equal = TRUE, alpha = 0.05){

  if (choose.model) {
    model.equal <- lme(y~met-1, data = data.long, random = ~1|id, method = "ML")
    model.unequal <- try(lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                             method = "ML"), silent = TRUE)
    # try opt = 'optim'
    if (inherits(model.unequal, "try-error")) {
      model.unequal <- try(lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                               method = "ML", control = lmeControl(opt = 'optim')), silent = TRUE)
      warning("Choudhary's parametric model with unequal variances fitted with lmeControl argument opt = 'optim'")
    }
    # try method = 'REML'
    if (inherits(model.unequal, "try-error")) {
      model.unequal <- try(lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                           method = "REML"), silent = TRUE)
      warning("Choudhary's parametric model with unequal variances fitted with REML method instead of ML")
    }

    # error if AIC cannot be compared
    if (inherits(try(AIC(model.unequal), silent = TRUE), "try-error")) {
      stop("Not able to compare AIC of Choudhary's parametric model with unequal variances. Consider rerunning with choose.model.ch.p = FALSE and var.equal = TRUE")
    }

    if (AIC(model.equal) <= AIC(model.unequal)) {
      model <- model.equal
      var.equal <- TRUE
    } else{
      model <- model.unequal
      var.equal <- FALSE
    }

    # error if model chose with unequal variances but returns an error
    if (!var.equal) {
      if (inherits(model.unequal, "try-error")) {
        stop("Not able to fit Choudhary's parametric model with unequal variances. Consider rerunning with choose.model.ch.p = FALSE and var.equal = TRUE")
      }
    }
  } else{
    if (var.equal) {
      model <- lme(y~met-1, data = data.long, random = ~1|id, method = "ML")
    } else{
      model <- try(lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                       method = "ML"), silent = TRUE)
      # try opt = 'optim'
      if (inherits(model, "try-error")) {
        model <- try(lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                         method = "ML", control = lmeControl(opt = 'optim')), silent = TRUE)
        warning("Choudhary's parametric model with unequal variances fitted with lmeControl argument opt = 'optim'")
      }
      # try method = 'REML'
      if (inherits(model, "try-error")) {
        model <- try(lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                         method = "REML"), silent = TRUE)
        warning("Choudhary's parametric model with unequal variances fitted with REML method instead of ML")
      }
      # error
      if (inherits(model, "try-error")) {
        stop("Not able to fit Choudhary's parametric model with unequal variances. Consider rerunning with var.equal = TRUE")
      }
    }
  }

  if (var.equal) {
    mu1 <- model$coefficients$fixed[1]
    mu2 <- model$coefficients$fixed[2]
    sigma <- model$sigma

    # try minAbsParApVar = 0.1
    if (is.na(model$apVar[4])) {
      model <- lme(y~met-1, data = data.long, random = ~1|id, method = "ML", control = lmeControl(minAbsParApVar = 0.1))
      warning("Choudhary's parametric model with equal variances fitted with lmeControl argument minAbsParApVar = 0.1")
    }
    if (is.na(model$apVar[4])) {
      stop("Choudhary's parametric model with equal variances did not converge")
    }

    # invI
    invI <- matrix(c(model$varFix[1], model$varFix[2], 0,
                     model$varFix[3], model$varFix[4], 0,
                     0, 0, model$apVar[4]),
                   nrow = 3, ncol = 3)

    tdi.hat <- TDI.Choudhary.P.EqVar(p, mu1, mu2, sigma, invI, ub, alpha)
  } else{
    mu1 <- model$coefficients$fixed[1]
    mu2 <- model$coefficients$fixed[2]
    sigma.1 <- model$sigma
    sigma.2 <- sigma.1*exp(as.matrix(model$modelStruct$varStruct))

    # try minAbsParApVar = 0.1
    if (is.character(model$apVar)) {
      model <- lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                   method = "ML", control = lmeControl(minAbsParApVar = 0.1))
      warning("Choudhary's parametric model with unequal variances fitted with lmeControl argument minAbsParApVar = 0.1")
    }
    if (is.character(model$apVar)) {
      stop("Choudhary's parametric model with unequal variances did not converge")
    }

    # invI
    zeros <- matrix(rep(0,4), nrow = 2)
    part1invI <- rbind(model$varFix, zeros)
    part2invI <- rbind(zeros, matrix(c(model$apVar[3, 3], model$apVar[2, 3] + model$apVar[3, 3],
                                       model$apVar[2, 3] + model$apVar[3, 3], sum(model$apVar[2:3, 2:3])),
                                     nrow = 2))
    invI <- cbind(part1invI, part2invI)
    colnames(invI) <- rownames(invI) <- c("mu1", "mu2", "log sigma1", "log sigma2")
    tdi.hat <- TDI.Choudhary.P.UneqVar(p, mu1, mu2, sigma.1, sigma.2, invI, ub, alpha)
  }

  out <- list(result = tdi.hat, fitted.model = model)
  return(out)
}
