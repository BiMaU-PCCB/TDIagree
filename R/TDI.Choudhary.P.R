#' TDI.Choudhary.P
#'
#' @param data.long name of the dataset containing at least 3 columns (measurement, subject effect, rater effect).
#' @param y measurement column name.
#' @param id subject effect column name. Must be a factor.
#' @param met rater effect column name. Must be a factor.
#' @param rep replication effect column name. Must be a factor. When there are no replicas rep = NA
#'            Default is NA.
#' @param p A value or vector of the percentiles for the TDI.
#' @param choose.model logical asking if the minimum AIC model should be chosen no matter the value of `var.equal`.
#' @param var.equal logical asking if there is homoscedasticity between methods.
#'                  Default is TRUE.
#' @param alpha confidence level for the TDI's upper bounds.
#'              Default is 0.05.
#'
#' @importFrom stats qchisq dnorm qnorm AIC
#' @importFrom nlme lme varIdent lmeControl
#'
#' @noRd


# auxiliar functions for TDI.Choudhary.P

TDI.Choudhary.P.EqVar <- function(p, mu1, mu2, sigma, invI, alpha = 0.05){
  mu.hat <- mu1 - mu2
  sigma.hat <- sqrt(2)*sigma
  q.hat <- sigma.hat*sqrt(qchisq(p, df = 1, ncp = (mu.hat/sigma.hat)^2))
  z.l <- -(q.hat + mu.hat)/sigma.hat
  z.u <- (q.hat - mu.hat)/sigma.hat
  s <- dnorm(z.l) + dnorm(z.u)
  partial.mu1 <- (1/(s*q.hat))*(dnorm(z.u) - dnorm(z.l))
  partial.mu2 <- (1/(s*q.hat))*(dnorm(z.l) - dnorm(z.u))
  partial.logsigma <- sigma.hat^2/(s*q.hat*sigma.hat)*(z.u*dnorm(z.u)-z.l*dnorm(z.l))
  G <- c(partial.mu1, partial.mu2, partial.logsigma)
  ul.TDI <- exp(log(q.hat) + qnorm(1-alpha)*sqrt(t(G)%*%as.matrix(invI)%*%G))
  res <- c(q.hat, ul.TDI)
  return(res)
}

TDI.Choudhary.P.UneqVar <- function(p, mu1, mu2, sigma1, sigma2, invI, alpha = 0.05){
  mu.hat <- mu1 - mu2
  sigma.hat <- sqrt(sigma1^2 + sigma2^2)
  q.hat <- sigma.hat*sqrt(qchisq(p, df = 1, ncp = (mu.hat/sigma.hat)^2))
  z.l <- -(q.hat + mu.hat)/sigma.hat
  z.u <- (q.hat - mu.hat)/sigma.hat
  s <- dnorm(z.l) + dnorm(z.u)
  partial.mu1 <- (1/(s*q.hat))*(dnorm(z.u) - dnorm(z.l))
  partial.mu2 <- (1/(s*q.hat))*(dnorm(z.l) - dnorm(z.u))
  partial.logsigma1 <- sigma1^2/(s*q.hat*sigma.hat)*(z.u*dnorm(z.u)-z.l*dnorm(z.l))
  partial.logsigma2 <- sigma2^2/(s*q.hat*sigma.hat)*(z.u*dnorm(z.u)-z.l*dnorm(z.l))
  G <- c(partial.mu1, partial.mu2, partial.logsigma1, partial.logsigma2)
  ul.TDI <- exp(log(q.hat) + qnorm(1-alpha)*sqrt(t(G)%*%as.matrix(invI)%*%G))
  res <- c(q.hat, ul.TDI)
  return(res)
}



TDI.Choudhary.P <- function(data.long, y, id, met, rep = NA, p,
                            choose.model = T, var.equal = T, alpha = 0.05){

  if(choose.model){
    model.equal <- lme(y~met-1, data = data.long, random = ~1|id, method = "ML")
    model.unequal <- try(lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                             method = "ML", control = lmeControl(opt = 'optim')), silent = TRUE)
    if(inherits(model.unequal, "try-error")){
      model.unequal <- lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                           method = "REML", control = lmeControl(opt = 'optim'))
      warning("Choudhary's parametric model with unequal variances fitted with REML method instead of ML")
    }
    if(AIC(model.equal) <= AIC(model.unequal)){
      model <- model.equal
      var.equal <- TRUE
    } else{
      model <- model.unequal
      var.equal <- FALSE
    }
  } else{
    if(var.equal){
      model <- lme(y~met-1, data = data.long, random = ~1|id, method = "ML")
    } else{
      model <- try(lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                       method = "ML", control = lmeControl(opt = 'optim')), silent = TRUE)
      if(inherits(model, "try-error")){
        model <- lme(y~met-1, data = data.long, random = ~1|id, weights = varIdent(form = ~1|met),
                     method = "REML", control = lmeControl(opt = 'optim'))
        warning("Choudhary's parametric model with unequal variances fitted with REML method instead of ML")
      }
    }
  }

  if (var.equal){
    mu1 <- model$coefficients$fixed[1]
    mu2 <- model$coefficients$fixed[2]
    sigma <- model$sigma

    # invI
    if(is.na(model$apVar[4])){
      # if the model does not converge, Var(log(sigma)) = 1/2*(number of subjects-2)
      invI <- matrix(c(model$varFix[1], model$varFix[2], 0,
                       model$varFix[3], model$varFix[4], 0,
                       0, 0, 1/(2*(model$dims$N-2))),
                     nrow = 3, ncol = 3)
    } else{
      invI <- matrix(c(model$varFix[1], model$varFix[2], 0,
                       model$varFix[3], model$varFix[4], 0,
                       0, 0, model$apVar[4]),
                     nrow = 3, ncol = 3)
    }
    TDI <- TDI.Choudhary.P.EqVar(p, mu1, mu2, sigma, invI, alpha)
  } else{
    mu1 <- model$coefficients$fixed[1]
    mu2 <- model$coefficients$fixed[2]
    sigma.1 <- model$sigma
    sigma.2 <- sigma.1*exp(as.matrix(model$modelStruct$varStruct))

    # invI
    zeros <- matrix(rep(0,4), nrow = 2)
    part1invI <- rbind(model$varFix, zeros)
    part2invI <- rbind(zeros, matrix(c(model$apVar[3, 3], model$apVar[2, 3] + model$apVar[3, 3],
                                       model$apVar[2, 3] + model$apVar[3, 3], sum(model$apVar[2:3, 2:3])),
                                     nrow = 2))
    invI <- cbind(part1invI, part2invI)
    colnames(invI) <- rownames(invI) <- c("mu1", "mu2", "log sigma1", "log sigma2")
    TDI <- TDI.Choudhary.P.UneqVar(p, mu1, mu2, sigma.1, sigma.2, invI, alpha)
  }

  out <- list(result = TDI, fitted.model = model)
  return(out)
}
