#' TDI estimation and inference
#'
#' @description
#' This function implements the estimation of the TDI and its corresponding \eqn{100(1-\alpha)\%} upper bound (UB),
#' where \eqn{\alpha} is the significance level, using the methods proposed by Choudhary (2008),
#' Escaramis \emph{et al.} (2010), Choudhary (2010) and Perez-Jaume and Carrasco (2015) in the case of two raters.
#' See \strong{Details} and \strong{References} for further information about these methods.
#'
#' @details
#' The methods of Choudhary (2008) and Escaramis \emph{et al.} (2010) are parametric methods based on linear mixed models
#' that assume normality of the data and linearity between the response and the effects (subjects, raters and
#' random errors). The methods of Choudhary (2010) and Perez-Jaume and Carrasco (2015) are non-parametric methods
#' based on the estimation of quantiles of the absolute value of the differences between raters. Non-parametric
#' methods are recommended when dealing with skewed data or other non-normally distributed data, such as count data.
#' In situations of normality, parametric methods are recommended. See \strong{References} for further details.
#'
#' @usage TDI(data, y, id, met, rep = NA,
#'     method = c("Choudhary P", "Escaramis", "Choudhary NP", "Perez-Jaume and Carrasco"),
#'     p = 0.9, dec.p = 2, dec.est = 3, tol = 10^(-8),
#'     choose.model.ch.p = TRUE, var.equal = TRUE, choose.model.es = TRUE, int = FALSE,
#'     alpha = 0.05, type = 8, R = 10000)
#'
#' @param data name of the dataset containing at least 3 columns (measurement, subject effect, rater effect).
#' @param y measurement column name.
#' @param id subject effect column name. Must be a factor.
#' @param met rater effect column name. Must be a factor.
#' @param rep replicate effect column name. Must be a factor. When there are no replicates use \code{rep = NA}. \cr
#'            Default is NA.
#' @param method name of the method(s) to estimate the TDI and UB. Options: \code{"Lin"} (Lin, 2000),
#'               \code{"Choudhary P"} (Choudhary, 2008), \code{"Escaramis"} (Escaramis \emph{et al.}, 2010),
#'               \code{"Choudhary NP"} (Choudhary, 2010) and \code{"Perez-Jaume and Carrasco"} (Perez-Jaume and Carrasco, 2015). \cr
#'               Default is \code{c("Choudhary P", "Escaramis", "Choudhary NP", "Perez-Jaume and Carrasco")}.
#' @param p a value or vector of the probabilities for estimation of the TDI, where \eqn{0<p<1}. Commonly, \eqn{p\geq 0.80}. \cr
#'          Default is 0.90.
#' @param dec.p number of decimals to display for \code{p} in the method \code{\link[TDIagree]{print.tdi}}. \cr
#'              Default is 2.
#' @param dec.est number of decimals to display for the estimates values in the method \code{\link[TDIagree]{print.tdi}}.
#'                Now implemented to display up to 4 decimals. In future versions, users will be able to choose to display
#'                more decimal places. \cr
#'                Default is 3.
#' @param tol tolerance to be used in the method of Escaramis \emph{et al.} (2010). \cr
#'            Default is 10^(-8).
#' @param choose.model.ch.p in the parametric method of Choudhary (2008) two methods can be fit, one with equal residual homoscedasticity
#'                          between raters and one with unequal residual homoscedasticity. This argument, if \code{TRUE}, chooses the
#'                          model with the minimum AIC. If \code{FALSE}, the argument \code{var.equal} must be specified. \cr
#'                          Default is \code{TRUE}.
#' @param var.equal logical asking if there is residual homoscedasticity between raters to choose the model in the parametric method of Choudhary (2008).
#'                  Ignored if \code{choose.model.ch.p} is \code{TRUE}. \cr
#'                  Default is \code{TRUE}.
#' @param choose.model.es in the method of Escaramis \emph{et al.} (2010) two methods can be fit, one including the subject--rater interaction
#'                        and one that does not. The model with interaction only applies to data with replicates. This argument, if \code{TRUE}, chooses the
#'                        model with the minimum AIC. If \code{FALSE}, the argument \code{int} must be specified. \cr
#'                        Default is \code{TRUE}.
#' @param int logical asking if there is interaction between subjects and methods to choose the model in the method of Escaramis \emph{et al.} (2010).
#'            The model with interaction only applies to data with replicates. Ignored if \code{choose.model.es} is \code{TRUE}. \cr
#'            Default is \code{FALSE}.
#' @param alpha confidence level for inference on the TDI. \cr
#'              Default is 0.05.
#' @param type in the method of Perez-Jaume and Carrasco (2015), a quantile is calculated to obtain the estimation of the tDI. This argument is an integer
#'             between 1 and 9 selecting one of the nine quantile algorithms (to be passed to \code{\link[stats]{quantile}}).
#'             Recommended 8 for continuous data and 3 for discrete data. \cr
#'             Default is 8.
#' @param R in the method of Perez-Jaume and Carrasco (2015), bootstrap is used for the estimation of the UB.
#'          This argument chooses the number of bootstrap replicates (to be passed to \code{\link[boot]{boot}}). \cr
#'          Default is 10000.
#'
#' @returns An object of class \code{tdi}, which is a list with five components:
#' \describe{
#'   \item{\code{result}}{an object of class \code{data.frame} with the TDI estimates and UBs of the methods specified for every probability.}
#'   \item{\code{fitted.models}}{a list with the fitted models of the parametric methods of Choudhary (2008) and Escaramis \emph{et al.} (2010).}
#'   \item{\code{params}}{a list with the values \code{dec.est} and \code{dec.p} to be used in the method \code{\link[TDIagree]{print.tdi}}.}
#'   \item{\code{data.long}}{an object of class \code{data.frame} with columns y, id, met (and rep if it applies) with the values of the measurement, subject identifiers,
#'                           rater (and replicate number if it applies) from the original data frame provided.}
#'   \item{\code{data.wide}}{an object of class \code{data.frame} with either:
#'                           \itemize{
#'                            \item{columns id, y_met1, y_met2 (in the case of no replicates) with the measurements of each method.}
#'                            \item{columns id, y_met1rep1, ..., y_met1rep\eqn{m}, y_met2rep1, ..., y_met2rep\eqn{m}, with the measurements of each method and each replicate, where \eqn{m} is the number of replicates.}}
#'                           Numbers 1 and 2 after met correspond to the first and second level of the column met in data, respectively.
#'                           Numbers 1, ..., \eqn{m} after rep correspond to the first, ..., \eqn{m}-th level of the column rep in data, respectively.}
#' }
#'
#'
#' @examples
#' set.seed(2025)
#'
#' n <- 100
#' y_A1 <- rnorm(n, 50, 10) # rater A, replicate 1
#' y_A2 <- rnorm(n, 50, 10) # rater A, replicate 2
#' y_B1 <- rnorm(n, 30, 15) # rater B, replicate 1
#' y_B2 <- rnorm(n, 30, 15) # rater B, replicate 2
#'
#' ex_data <- data.frame(y = c(y_A1, y_A2, y_B1, y_B2), rater = factor(rep(c("A", "B"), each = 2*n)),
#'                       replicate = factor(rep(rep(1:2, each = n), 2)), subj = factor(rep(1:n, 4)))
#'
#' tdi <- TDI(ex_data, y, subj, rater, replicate, p = c(0.8, 0.9),
#'            method = c("Choudhary P", "Escaramis", "Choudhary NP", "Perez-Jaume and Carrasco"))
#' tdi$result
#' tdi$fitted.models
#' tdi$data.long
#' tdi$data.wide
#'
#' @seealso \code{\link[TDIagree]{print.tdi}}, \code{\link[TDIagree]{plot.tdi}}
#'
#' @references Efron, B., & Tibshirani, R. (1993). An Introduction to the Bootstrap; Chapman and 913 Hall. Inc.: New York, NY, USA, 914.
#'
#'    Lin, L. I. K. (2000). Total deviation index for measuring individual agreement with applications in laboratory performance and bioequivalence. Statistics in medicine, 19.2:255-270.
#'
#'    Choudhary, P. K. (2008). A tolerance interval approach for assessment of agreement in method comparison studies with repeated measurements. Journal of Statistical Planning and inference, 138(4), 1102-1115.
#'
#'    Escaramís, G., Ascaso, C., & Carrasco, J. L. (2010). The total deviation index estimated by tolerance intervals to evaluate the concordance of measurement devices. BMC Medical Research Methodology, 10, 1-12.
#'
#'    Choudhary, P. K. (2010). A unified approach for nonparametric evaluation of agreement in method comparison studies. The International Journal of Biostatistics, 6(1).
#'
#'    Perez‐Jaume, S., & Carrasco, J. L. (2015). A non‐parametric approach to estimate the total deviation index for non‐normal data. Statistics in medicine, 34(25), 3318-3335.
#'
#' @export

TDI <- function(data, y, id, met, rep = NA,
                method = c("Choudhary P", "Escaramis", "Choudhary NP", "Perez-Jaume and Carrasco"),
                p = 0.9, dec.p = 2, dec.est = 3, tol = 10^(-8),
                choose.model.ch.p = TRUE, var.equal = TRUE, choose.model.es = TRUE, int = FALSE,
                alpha = 0.05, type = 8, R = 10000) {

  # columns
  y <- paste(deparse(substitute(y)), collapse = "")
  y <- gsub("  ", "", y)
  id <- paste(deparse(substitute(id)), collapse = "")
  id <- gsub("  ", "", id)
  met <- paste(deparse(substitute(met)), collapse = "")
  met <- gsub("  ", "", met)
  rep <- paste(deparse(substitute(rep)), collapse = "")
  rep <- gsub("  ", "", rep)

  # stop checks
  if (!is.data.frame(data)){
    stop("data must be a data.frame")
  }

  if (rep != "NA") {
    if (!all(c(y, id, met, rep) %in% names(data))){
      stop("'y', 'id', 'met' and 'rep' must be columns in 'data'")
    }
    m =  length(levels(data[, rep]))
    if (any(duplicated(c(y, id, met, rep)))){
      stop("two of the column identifiers are the same")
    }
    if (!is.factor(data[, id]) | !is.factor(data[, met])|!is.factor(data[, rep])){
      stop("'id', 'met' and 'rep' columns must be factors")
    }
  }
  else{
    rep <- NA
    if (!all(c(y, id, met) %in% names(data))){
      stop("'y', 'id' and 'met' must be columns in 'data'")
    }
    m = 1
    if (any(duplicated(c(y, id, met)))){
      stop("two of the column identifiers are the same")
    }
    if (!is.factor(data[, id]) | !is.factor(data[, met])){
      stop("'id' and 'met' columns must be factors")
    }
  }

  if ((var.equal != FALSE) & (var.equal != TRUE)){
    stop("'var.equal' must be FALSE (model with unequal variances) or TRUE (model with equal variances)")
  }

  if ((int != FALSE) & (int != TRUE)){
    stop("'int' must be FALSE (model with no interaction) or TRUE (model with interaction)")
  }

  if (!all(table(data[, met]) == table(data[, met])[1])) {
    stop("the design must be balanced, there must be the same number of measurements per method")
  }

  # warning checks
  if (!is.na(rep)){
    if (length(levels(data[, rep])) == 1){
      warning("'rep' column has only 1 level. Consider omitting 'rep' call. Calculation proceeds assuming rep = NA")
      rep <- NA
    }
  }

  if (length(levels(data[, id])) < 6){
    warning("number of subjects less than 6. Proceed with caution")
  }

  if ((int == TRUE) & is.na(rep)){
    warning("no replicates specified in 'rep', unable to fit model with interaction. Calculation proceeds assuming int = FALSE")
    int <- FALSE
  }

  # extracting necessary columns of data
  if (!is.na(rep)) {
    data.long <- data.frame(y = data[, y], id = data[, id],
                        met = data[, met], rep = data[, rep])
    data.long$z <- interaction(data.long$met, data.long$rep)
    data.wide <- data.frame(id = unique(data.long$id))
    for (j in 1:(length(levels(data.long$z))/2)){
      data.wide[, paste0("y_met1rep", j)] <- data.long[data.long$z == levels(data.long$z)[2*j-1], "y"]
      data.wide[, paste0("y_met2rep", j)] <- data.long[data.long$z == levels(data.long$z)[2*j], "y"]
    }
  } else{
    data.long <- data.frame(y = data[, y], id = data[, id], met = data[, met])
    y1 <- subset(data.long, met==levels(data.long[, "met"])[1])$y
    y2 <- subset(data.long, met==levels(data.long[, "met"])[2])$y
    data.wide <- data.frame(id=unique(data.long$id), y_met1 = y1, y_met2 = y2)
  }

  # defining outputs
  result <- data.frame(p = p,
                       tdi.lin.p = NA, ub.lin.p = NA,
                       tdi.ch.p = NA, ub.ch.p = NA,
                       tdi.es.p = NA, ub.es.p = NA,
                       tdi.ch.np = NA, ub.ch.np = NA,
                       tdi.pc.np = NA, ub_p.pc.np = NA, ub_nlog.pc.np = NA)

  output <- list()

  model.ch <- "Choudhary's parametric method not called"
  model.es <- "Escaramis' method not called"

  for (i in 1:length(p)){
    # method of Lin
    if ("Lin" %in% method) {
      warning("Lin's method to be implemented in future versions")
    }

    # parametric method of Choudhary
    if ("Choudhary P" %in% method) {
      choud.p <- TDI.Choudhary.P(data.long, y, id, met, rep, p = p[i],
                                 choose.model = choose.model.ch.p, var.equal = var.equal, alpha = alpha)
      result$tdi.ch.p[i] <- choud.p$result[1]
      result$ub.ch.p[i] <- choud.p$result[2]
      model.ch <- choud.p$fitted.model
    }

    # method of Escaramis et al.
    if ("Escaramis" %in% method) {
      escaramis <- TDI.Escaramis(data.long, y, id, met, rep, p = p[i],
                                 choose.model = choose.model.es, tol = tol, int = int, alpha = alpha)
      result$tdi.es.p[i] <- escaramis$result[1]
      result$ub.es.p[i] <- escaramis$result[2]
      model.es <- escaramis$fitted.model
    }

    # non-parametric method of Choudhary
    if ("Choudhary NP" %in% method) {
      choud.np <- ChoudNonPar(data, y, id, met, rep, p = p[i], alpha = alpha)
      result$tdi.ch.np[i] <- choud.np[1]
      result$ub.ch.np[i] <- choud.np[2]
    }

    # method of Perez-Jaume and Carrasco
    if ("Perez-Jaume and Carrasco" %in% method) {
      if (!is.na(rep)){
        for (k in 1:length(levels(data.long$rep))) {
          for (l in 1:length(levels(data.long$rep))) {
            dif <- data.wide[, paste0("y_met1rep", k)] - data.wide[, paste0("y_met2rep", l)]
            if (k == 1 & l == 1) {
              d <- dif
            } else{
              d <- append(d, dif)
            }
          }
        }
      } else{
        d <- data.wide$y_met1 - data.wide$y_met2
      }
      nonpar <- NonParTDI(d, p[i], alpha, type, R)
      result$tdi.pc.np[i] <- nonpar[1]
      result$ub_p.pc.np[i] <- nonpar[2]
      result$ub_nlog.pc.np[i] <- nonpar[3]
    }
  }

  result <- result[ , colSums(is.na(result))==0]

  output$result <- result
  output$fitted.models$choudhary <- model.ch
  output$fitted.models$escaramis <- model.es
  output$params$dec.est <- dec.est
  output$params$dec.p <- dec.p
  if (is.na(rep)){
    output$data.long <- data.long
  } else{
    output$data.long <- data.long[, 1:4]
  }

  output$data.wide <- data.wide

  class(output) <- "tdi"

  return(output)
}
