#' TDI estimation and inference
#'
#' @description
#' This function implements the estimation of the TDI and its corresponding \eqn{100(1-\alpha)\%} upper bound (UB),
#' where \eqn{\alpha} is the significance level, using the methods proposed by Choudhary (2007),
#' Escaramis \emph{et al.} (2010), Choudhary (2010) and Perez-Jaume and Carrasco (2015) in the case of two raters.
#' See \strong{Details} and \strong{References} for further information about these methods.
#'
#' @details
#' The methods of Choudhary (2007) and Escaramis \emph{et al.} (2010) are parametric methods based on linear mixed models
#' that assume normality of the data and linearity between the response and the effects (subjects, raters and
#' random errors). The linear mixed models are fitted using the function \code{\link[nlme]{lme}} from the \code{nlme} package.
#' The methods of Choudhary (2010) and Perez-Jaume and Carrasco (2015) are non-parametric methods
#' based on the estimation of quantiles of the absolute value of the differences between raters. Non-parametric
#' methods are recommended when dealing with skewed data or other non-normally distributed data, such as count data.
#' In situations of normality, parametric methods are recommended. See \strong{References} for further details.
#'
#' @usage TDI(data, y, id, met, rep = NA,
#'     method = c("Choudhary P", "Escaramis et al.",
#'                "Choudhary NP", "Perez-Jaume and Carrasco"),
#'     p = 0.9, ub = TRUE, boot.type = c("differences", "cluster"),
#'     type = 8, R = 10000, dec.p = 2, dec.est = 3,
#'     choose.model.ch.p = TRUE, var.equal = TRUE,
#'     choose.model.es = TRUE, int = FALSE, tol = 10^(-8), add.es = NULL,
#'     alpha = 0.05)
#'
#' @param data name of the dataset, of class \code{data.frame}, containing at least 3 columns (quantitative measurement, subject effect, rater effect).
#' @param y quantitative measurement column name.
#' @param id subject effect column name. The corresponding column of \code{data} must be a factor.
#' @param met rater effect column name. The corresponding column of \code{data} must be a factor.
#' @param rep replicate effect column name. When there are no replicates the user should use \code{rep = NA}. When there are replicates, the corresponding column of data must be a factor. \cr
#'            The default value is \code{NA}.
#' @param method name of the method(s) to estimate the TDI and UB. The options are:
#'               \code{"Choudhary P"} (Choudhary, 2007), \code{"Escaramis et al."} (Escaramis \emph{et al.}, 2010),
#'               \code{"Choudhary NP"} (Choudhary, 2010) and \code{"Perez-Jaume and Carrasco"} (Perez-Jaume and Carrasco, 2015).
#'               This argument is not case-sensitive and is passed to \code{\link[base]{match.arg}}. \cr
#'               The default value is \code{c("Choudhary P", "Escaramis et al.", "Choudhary NP", "Perez-Jaume and Carrasco")}, so all approaches are executed by default.
#' @param p a value or vector of the proportion(s) for estimation of the TDI, where \eqn{0<p<1}. Commonly, \eqn{p\geq 0.80}. \cr
#'          The default value is 0.90.
#' @param ub logical asking whether the UBs should be computed. \cr
#'           The default value is \code{TRUE}.
#' @param boot.type name of the bootstrap approach(es) to be used in the method of Perez-Jaume and Carrasco (2015). There are two
#'                  different options when there are replicates: to bootstrap the vector of the within-subject differences (\code{"differences"})
#'                  or to bootstrap at subject level (\code{"cluster"}). This is, not all the differences coming from the same subject
#'                  need to be bootstrapped together in the first one but all the measurements from the same subjects have to be bootstrapped
#'                  together in the second one. This argument is passed to \code{\link[base]{match.arg}} \cr
#'                  The default value is \code{c("differences", "cluster")}, so all approaches are executed by default.
#' @param type in the method of Perez-Jaume and Carrasco (2015), a quantile is calculated to obtain the estimation of the TDI. This argument is an integer
#'             between 1 and 9 selecting one of the nine quantile algorithms (to be passed to \code{\link[stats]{quantile}}).
#'             We recommend 8 for continuous data and 3 for discrete data. \cr
#'             The default value is 8.
#' @param R in the method of Perez-Jaume and Carrasco (2015), bootstrap is used for the estimation of the UB.
#'          This argument chooses the number of bootstrap replicates (to be passed to \code{\link[boot]{boot}}). \cr
#'          The default value is 10000.
#' @param dec.p number of decimals to display for \code{p} in the method \code{\link[TDIagree]{print.tdi}}. \cr
#'              The default value is 2.
#' @param dec.est number of decimals to display for the estimates in the method \code{\link[TDIagree]{print.tdi}}. Up to 4 decimals. \cr
#'                The default value is 3.
#' @param choose.model.ch.p in the parametric method of Choudhary (2007) two methods can be fit, one with equal residual homoscedasticity
#'                          between raters and one with unequal residual homoscedasticity. This argument, if \code{TRUE}, chooses the
#'                          model with the minimum AIC. If \code{FALSE}, the argument \code{var.equal} must be specified. \cr
#'                          The default value is \code{TRUE}.
#' @param var.equal logical asking if there is residual homoscedasticity between raters to choose the model in the parametric method of Choudhary (2007).
#'                  If \code{choose.model.ch.p} is set to \code{TRUE}, this argument is ignored. \cr
#'                  The default value is \code{TRUE}.
#' @param choose.model.es in the method of Escaramis \emph{et al.} (2010) two methods can be fit, one including the subject--rater interaction
#'                        and one that does not. The model with interaction only applies to data with replicates. This argument, if \code{TRUE}, chooses the
#'                        model with the minimum AIC. If \code{FALSE}, the argument \code{int} must be specified. \cr
#'                        The default value is \code{TRUE}.
#' @param int logical asking if there is interaction between subjects and methods to choose the model in the method of Escaramis \emph{et al.} (2010).
#'            The model with interaction only applies to data with replicates. If \code{choose.model.es} is set to \code{TRUE}, this argument is ignored. \cr
#'            The default value is \code{FALSE}.
#' @param tol tolerance to be used in the method of Escaramis \emph{et al.} (2010). \cr
#'            The default value is 10^(-8).
#' @param add.es name of the columns in \code{data} that will be added to the model (as fixed effects) of the method of Escaramis \emph{et al.} (2010). It must be passed as a column
#'               name or vector of column names. \cr
#'               The default value, \code{NULL}, indicates that no extra variables are to be added in the model.
#' @param alpha significance level for inference on the TDI. \cr
#'              The default value is 0.05.
#'
#' @returns An object of class \code{tdi}, which is a list with five components:
#' \describe{
#'   \item{\code{result}}{an object of class \code{data.frame} with the TDI estimates and UBs of the methods specified for every proportion.}
#'   \item{\code{fitted.models}}{a list with the fitted models of the parametric methods of Choudhary (2007) and Escaramis \emph{et al.} (2010).}
#'   \item{\code{params}}{a list with the values \code{dec.est}, \code{dec.p}, \code{ub}, \code{method} and \code{alpha} to be used in the method \code{\link[TDIagree]{print.tdi}}
#'                        and in the method \code{\link[TDIagree]{plot.tdi}}.}
#'   \item{\code{data.long}}{an object of class \code{data.frame} with columns y, id, met (and rep if it applies) with the values of the measurement, subject identifiers,
#'                           rater (and replicate number if it applies) from the original data frame provided.}
#'   \item{\code{data.wide}}{an object of class \code{data.frame} with either:
#'                           \itemize{
#'                            \item{columns id, y_met1, y_met2 (in the case of no replicates) with the measurements of each method.}
#'                            \item{columns id, y_met1rep1,..., y_met1rep\eqn{m}, y_met2rep1,..., y_met2rep\eqn{m}, with the measurements of each method and each replicate, where \eqn{m} is the number of replicates.}}
#'                           Numbers 1 and 2 after met correspond to the first and second level of the column met in data, respectively.
#'                           Numbers 1,..., \eqn{m} after rep correspond to the first,..., \eqn{m}-th level of the column rep in data, respectively.}
#' }
#'
#'
#' @examples
#' # normal data, parametric methods more suitable
#'
#' set.seed(2025)
#'
#' n <- 100
#'
#' mu.ind <- rnorm(n, 0, 7)
#'
#' epsA1 <- rnorm(n, 0, 3)
#' epsA2 <- rnorm(n, 0, 3)
#' epsB1 <- rnorm(n, 0, 3)
#' epsB2 <- rnorm(n, 0, 3)
#'
#' y_A1 <- 50 + mu.ind + epsA1 # rater A, replicate 1
#' y_A2 <- 50 + mu.ind + epsA2 # rater A, replicate 2
#' y_B1 <- 40 + mu.ind + epsB1 # rater B, replicate 1
#' y_B2 <- 40 + mu.ind + epsB2 # rater B, replicate 2
#'
#' ex_data <- data.frame(y = c(y_A1, y_A2, y_B1, y_B2),
#'                       rater = factor(rep(c("A", "B"), each = 2*n)),
#'                       replicate = factor(rep(rep(1:2, each = n), 2)),
#'                       subj = factor(rep(1:n, 4)))
#'
#' tdi <- TDI(ex_data, y, subj, rater, replicate, p = c(0.8, 0.9),
#'            method = c("Choudhary P", "Escaramis et al.",
#'                       "Choudhary NP", "Perez-Jaume and Carrasco"),
#'            boot.type = "cluster", R = 1000)
#'
#' tdi$result
#' tdi$fitted.models
#' tdi$data.long
#' tdi$data.wide
#'
#'
#' # non-normal data, non-parametric methods more suitable
#'
#' tdi.aml <- TDI(AMLad, mrd, id, met, rep, p = c(0.85, 0.95), boot.type = "cluster",
#'                dec.est = 4, R = 1000)
#' tdi.aml$result
#' tdi.aml$fitted.models
#' tdi.aml$data.long
#' tdi.aml$data.wide
#'
#' @seealso \code{\link[TDIagree]{print.tdi}}, \code{\link[TDIagree]{plot.tdi}}
#'
#' @references Efron, B., & Tibshirani, R. (1993). An Introduction to the Bootstrap; Chapman and Hall. Inc.: New York, NY, USA, 914.
#'
#'    Lin, L. I. K. (2000). Total deviation index for measuring individual agreement with applications in laboratory performance and bioequivalence. Statistics in Medicine, 19(2):255-270.
#'
#'    Choudhary, P. K. (2007). A tolerance interval approach for assessment of agreement with left censored data. Journal of Biopharmaceutical Statistics, 17(4), 583-594.
#'
#'    Escaramis, G., Ascaso, C., & Carrasco, J. L. (2010). The total deviation index estimated by tolerance intervals to evaluate the concordance of measurement devices. BMC Medical Research Methodology, 10, 1-12.
#'
#'    Choudhary, P. K. (2010). A unified approach for nonparametric evaluation of agreement in method comparison studies. The International Journal of Biostatistics, 6(1).
#'
#'    Perez‐Jaume, S., & Carrasco, J. L. (2015). A non‐parametric approach to estimate the total deviation index for non‐normal data. Statistics in Medicine, 34(25), 3318-3335.
#'
#' @export

TDI <- function(data, y, id, met, rep = NA,
                method = c("Choudhary P", "Escaramis et al.",
                           "Choudhary NP", "Perez-Jaume and Carrasco"),
                p = 0.9, ub = TRUE, boot.type = c("differences", "cluster"),
                type = 8, R = 10000, dec.p = 2, dec.est = 3,
                choose.model.ch.p = TRUE, var.equal = TRUE,
                choose.model.es = TRUE, int = FALSE, tol = 10^(-8), add.es = NULL,
                alpha = 0.05) {

  # columns
  y <- paste(deparse(substitute(y)), collapse = "")
  y <- gsub("  ", "", y)
  id <- paste(deparse(substitute(id)), collapse = "")
  id <- gsub("  ", "", id)
  met <- paste(deparse(substitute(met)), collapse = "")
  met <- gsub("  ", "", met)
  rep <- paste(deparse(substitute(rep)), collapse = "")
  rep <- gsub("  ", "", rep)
  add.es <- paste(deparse(substitute(add.es)), collapse = "")
  add.es <- gsub(" ", "", add.es)
  add.es <- unlist(strsplit(gsub("\\)", "", gsub("c\\(", "", add.es)), ","))
  if ("NULL" %in% add.es) {
    add.es <- NULL
  }

  # argument type checks
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }
  if (rep != "NA") {
    if (!all(c(y, id, met, rep) %in% names(data))) {
      stop("'y', 'id', 'met' and 'rep' must be columns in 'data'")
    }
    if (any(duplicated(c(y, id, met, rep)))) {
      stop("two of the column identifiers are the same")
    }
    if (!is.factor(data[, id]) | !is.factor(data[, met]) | !is.factor(data[, rep])) {
      stop("'id', 'met' and 'rep' columns must be factors")
    }
  } else{
    rep <- NA
    if (!all(c(y, id, met) %in% names(data))) {
      stop("'y', 'id' and 'met' must be columns in 'data'")
    }
    if (any(duplicated(c(y, id, met)))) {
      stop("two of the column identifiers are the same")
    }
    if (!is.factor(data[, id]) | !is.factor(data[, met])) {
      stop("'id' and 'met' columns must be factors")
    }
  }
  choices.method <- tolower(c("Choudhary P", "Escaramis et al.", "Choudhary NP", "Perez-Jaume and Carrasco"))
  method.pre <- tolower(method)
  method <- match.arg(method.pre, choices.method, several.ok = TRUE)
  if (length(method) != length(method.pre)) {
    warning("one of the methods called does not exist and will not be estimated")
  }
  p <- unique(p)
  if (!is.numeric(p)) {
    stop("'p' must be a numeric value or vector")
  }
  if (!is.logical(ub)) {
    stop("'ub' must be logical")
  }
  boot.type <- match.arg(boot.type, several.ok = TRUE)
  if (!is.numeric(type)) {
    stop("'type' must be numeric")
  } else if (type != 1 & type != 2 & type != 3 & type != 4 & type != 5 & type != 6 & type != 7 & type != 8 & type != 9) {
    stop("'type' must be an integer between 1 and 9")
  }
  if (!is.numeric(R)) {
    stop("'R' must be numeric")
  }
  if (!is.numeric(dec.p)) {
    stop("'dec.p' must be numeric")
  } else if (dec.p > 4) {
    stop("'dec.p' must be less or equal than 4")
  }
  if (!is.numeric(dec.est)) {
    stop("'dec.est' must be numeric")
  } else if (dec.est > 4) {
    stop("'dec.est' must be less or equal than 4")
  }
  if (!is.logical(choose.model.ch.p)) {
    stop("'choose.model.ch.p' must be logical")
  }
  if (!is.logical(var.equal)) {
    stop("'var.equal' must be logical")
  }
  if (!is.logical(choose.model.es)) {
    stop("'choose.model.es' must be logical")
  }
  if (!is.logical(int)) {
    stop("'int' must be logical")
  }
  if (!is.numeric(tol)) {
    stop("'tol' must be numeric")
  }
  if (!is.null(add.es)) {
    if (!all(add.es %in% names(data))) {
      stop("'add.es' must be NULL or columns in 'data'")
    }
  }
  if (!is.numeric(alpha)) {
    stop("'alpha' must be numeric")
  } else if (alpha <= 0 | alpha >= 1) {
    stop("'alpha' must be between 0 and 1, both not included")
  }

  # warning checks
  if (any(table(data[, id]) == 0)) {
    warning("'id' column has (at least) one empty level. Calculation proceeds ignoring the empty level(s)")
    data[, id] <- factor(data[, id])
  }
  if (any(table(data[, met]) == 0)) {
    warning("'met' column has (at least) one empty level. Calculation proceeds ignoring the empty level(s)")
    data[, met] <- factor(data[, met])
  }
  if (nlevels(data[, met]) != 2) {
    stop("this package only supports two levels in 'met' column")
  }
  if (!is.na(rep)) {
    if (any(table(data[, rep]) == 0)) {
      warning("'rep' column has (at least) one empty level. Calculation proceeds ignoring the empty level(s)")
      data[, rep] <- factor(data[, rep])
    }
    if (nlevels(data[, rep]) == 1) {
      warning("'rep' column only has one level. Consider omitting 'rep' call. Calculation proceeds assuming rep = NA")
      rep <- NA
    }
  }
  if (nlevels(data[, id]) < 6){
    warning("number of subjects less than 6. Proceed with caution")
  }
  if (int & is.na(rep)){
    warning("no replicates specified in 'rep', unable to fit model with interaction. Calculation proceeds assuming int = FALSE")
    int <- FALSE
  }
  if (ub) {
    if (is.na(rep) & "cluster" %in% boot.type){
      warning("no replicates specified in 'rep', cluster bootstrap and differences bootstrap are the same. Consider omitting 'cluster' in boot.type call. Calculation proceeds assuming boot.type = 'differences'")
      boot.type <- "differences"
    }
  }

  # stop checks
  if (!all(table(data[, met]) == table(data[, met])[1])) {
    stop("the design must be balanced, there must be the same number of measurements per method")
  }
  if (is.na(rep)){
    if (sum(table(data[, id], data[, met]) != 1) > 0) {
      stop("the design must be balanced, each subject should have one measurement per method")
    }
  } else{
    if (sum(table(data[, id], data[, met], data[, rep]) != 1) > 0) {
      stop("the design must be balanced, each subject should have one measurement per method and replicate")
    }
  }
  # extracting number of replicate
  if (!is.na(rep)) {
    m <- nlevels(data[, rep])
  } else{
    m <- 1
  }

  # extracting necessary columns of data
  if (!is.na(rep)) {
    data.long <- data.frame(y = data[, y], id = data[, id],
                            met = data[, met], rep = data[, rep])
    if (!is.null(add.es)) {
      for(i in 1:length(add.es)){
        data.long[, add.es[i]] <- data[, add.es[i]]
      }
    }
    data.long$z <- interaction(data.long$met, data.long$rep)
    data.wide <- data.frame(id = unique(data.long$id))
    for (j in 1:(nlevels(data.long$z)/2)){
      data.wide[, paste0("y_met1rep", j)] <- data.long[data.long$z == levels(data.long$z)[2*j-1], "y"]
      data.wide[, paste0("y_met2rep", j)] <- data.long[data.long$z == levels(data.long$z)[2*j], "y"]
    }
  } else{
    data.long <- data.frame(y = data[, y], id = data[, id], met = data[, met])
    if (!is.null(add.es)) {
      for(i in 1:length(add.es)){
        data.long[, add.es[i]] <- data[, add.es[i]]
      }
    }
    y1 <- subset(data.long, met==levels(data.long[, "met"])[1])$y
    y2 <- subset(data.long, met==levels(data.long[, "met"])[2])$y
    data.wide <- data.frame(id=unique(data.long$id), y_met1 = y1, y_met2 = y2)
  }

  # stop check NAs
  if (any(is.na(data.long))) {
    stop("NAs not supported for TDI estimation. Rerun without NAs")
  }

  # defining outputs
  result <- data.frame(p = p,
                       tdi.ch.p = NA, ub.ch.p = NA,
                       tdi.es.p = NA, ub.es.p = NA,
                       tdi.ch.np = NA, ub.ch.np = NA,
                       tdi.pc.np = NA,
                       ub_p_db.pc.np = NA, ub_n_db.pc.np = NA, ub_e_db.pc.np = NA, ub_b_db.pc.np = NA,
                       ub_p_cb.pc.np = NA, ub_n_cb.pc.np = NA, ub_e_cb.pc.np = NA, ub_b_cb.pc.np = NA)

  output <- list()

  model.ch <- "Parametric method of Choudhary not called"
  model.es <- "Method of Escaramis et al. not called"

  for (i in 1:length(p)) {

    # parametric method of Choudhary
    if ("choudhary p" %in% method) {
      choud.p <- try(TDI.Choudhary.P(data.long, p = p[i], ub = ub,
                                     choose.model = choose.model.ch.p, var.equal = var.equal, alpha = alpha), silent = FALSE)
      if (inherits(choud.p, "try-error") & length(method) == 1) {
        stop("Parametric method of Choudhary did not converge")
      } else if (inherits(choud.p, "try-error")) {
        warning("Parametric method of Choudhary did not converge")
        result$tdi.ch.p[i] <- NA
        result$ub.ch.p[i] <- NA
        model.ch <- "Parametric method of Choudhary did not converge"
      } else{
        result$tdi.ch.p[i] <- choud.p$result[1]
        result$ub.ch.p[i] <- choud.p$result[2]
        model.ch <- choud.p$fitted.model
      }
    }

    # method of Escaramis et al.
    if ("escaramis et al." %in% method) {
      escaramis <- TDI.Escaramis(data.long, rep, p = p[i], ub = ub,
                                 choose.model = choose.model.es, tol = tol, int = int,
                                 add.es = add.es, alpha = alpha)
      result$tdi.es.p[i] <- escaramis$result[1]
      result$ub.es.p[i] <- escaramis$result[2]
      model.es <- escaramis$fitted.model
    }

    # non-parametric method of Choudhary
    if ("choudhary np" %in% method) {
      choud.np <- TDI.Choudhary.NP(data.long, rep, p = p[i], ub = ub, alpha = alpha)
      result$tdi.ch.np[i] <- choud.np[1]
      result$ub.ch.np[i] <- choud.np[2]
    }

    # method of Perez-Jaume and Carrasco
    if ("perez-jaume and carrasco" %in% method) {
      pjc <- TDI.Perez_Jaume_and_Carrasco(data.wide, p = p[i], ub = ub, boot.type = boot.type,
                                          alpha = alpha, type = type, R = R)
      result$tdi.pc.np[i] <- pjc[1]
      result$ub_p_db.pc.np[i] <- pjc[2]
      result$ub_n_db.pc.np[i] <- pjc[3]
      result$ub_e_db.pc.np[i] <- pjc[4]
      result$ub_b_db.pc.np[i] <- pjc[5]
      result$ub_p_cb.pc.np[i] <- pjc[6]
      result$ub_n_cb.pc.np[i] <- pjc[7]
      result$ub_e_cb.pc.np[i] <- pjc[8]
      result$ub_b_cb.pc.np[i] <- pjc[9]
    }
  }

  result <- result[ , colSums(is.na(result))==0]

  output$result <- result
  output$fitted.models$choudhary <- model.ch
  output$fitted.models$escaramis <- model.es
  output$params$dec.est <- dec.est
  output$params$dec.p <- dec.p
  output$params$ub <- ub
  output$params$method <- method
  output$params$alpha <- alpha


  if (is.na(rep)) {
    output$data.long <- data.long[, 1:3]
  } else{
    output$data.long <- data.long[, 1:4]
  }

  output$data.wide <- data.wide

  class(output) <- "tdi"

  return(output)
}
