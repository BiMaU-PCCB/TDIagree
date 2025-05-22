#' Bland-Altman plot
#'
#' @description
#' This function creates a Bland-Altman plot from Altman and Bland (1983), which is used to evaluate the agreement among the quantitative measures taken by two raters.
#' The plot displays the mean of the measurements from both raters in the x-axis and the differences between the measures taken by the two raters
#' in the y-axis. It can also display the TDI and UB estimates from the call of the function \code{\link[TDIagree]{TDI}} as well as the limits of
#' agreement (LoA) from Bland and Altman (1986).
#'
#' @details
#' The LoA are computed using the formula \eqn{\bar{d}\pm z_{1-\frac{\alpha}{2}}\cdot \text{sd}(d)}, where \eqn{z_{1-\frac{\alpha}{2}}} is the \eqn{(1-\frac{\alpha}{2})}-th
#' quantile of the standard normal distribution, \eqn{d} is the vector containing the differences between the two raters and \eqn{\bar{d}} represents their mean.
#'
#' @param x input object of class \code{tdi} resulting from a call of the function \code{\link[TDIagree]{TDI}}.
#' @param tdi logical indicating whether the \eqn{\pm}TDI estimate(s) should be added to the plot as solid lines. \cr
#'            The default value is \code{FALSE}.
#' @param ub logical indicating whether the \eqn{\pm}UB estimate(s) should be added to the plot as dashed lines. \cr
#'           The default value is \code{FALSE}.
#' @param loa logical indicating whether the LoA should be added to the plot as dotted lines. \cr
#'            The default value is \code{FALSE}.
#' @param method name of the method(s) for which the TDI or the UB estimates will be added to the plot. If both \code{tdi} and \code{ub}
#'               are set to \code{FALSE}, this argument is ignored. This argument is not case-sensitive and is passed to \code{\link[base]{match.arg}}. \cr
#'               The default value, \code{NULL}, indicates that, for the measures specified, all the methods for which the TDI (and/or UB) has been
#'               computed in the call of the function \code{\link[TDIagree]{TDI}} are to be added to the plot.
#' @param ub.pc name of the technique for the estimated UB to be added from the method of Perez-Jaume and Carrasco (2015). Possible values are: \code{p_db}, \code{n_db}, \code{e_db}, \code{b_db},
#'              \code{p_cb}, \code{n_cb}, \code{e_cb} and \code{b_cb}. The bootstrap approach (differences or cluster) is indicated with "db" and "cb" and the
#'              strategy (based on percentiles, the normal distribution, the empirical method or the BC\eqn{_a}) is indicated with "p", "n", "e" and "b". \cr
#'              The default value, \code{NULL}, indicates that the first estimated UB is to be added to the plot.
#' @param p value of the proportion for which the TDI and/or UB (depending on the value of the arguments \code{tdi} and \code{ub})
#'          are to be added to the plot. If both \code{tdi} and \code{ub} are set to \code{FALSE}, this argument is ignored. \cr
#'          The default value, \code{NULL}, indicates that only the first proportion passed to the call of the function \code{\link[TDIagree]{TDI}}
#'          is to be considered.
#' @param loess logical indicating whether a smooth curve computed by \code{\link[stats]{loess.smooth}} should be added to the plot as a dotdashed curve. \cr
#'              The default value is \code{FALSE}.
#' @param method.col colour palette to be used in the drawing of TDIs and/or UBs. A colour should be indicated for every method asked. It is assumed that the colours
#'                   are passed in the same order as the methods passed to \code{method}. If both \code{tdi} and \code{ub} are set to \code{FALSE},
#'                   this argument is ignored. \cr
#'                   The default value, \code{NULL}, indicates that the following palette should be used:
#'                   \code{"#f3df6c"}, \code{"#9c964a"}, \code{"#f4b5bd"} and \code{"#85d4e3"} corresponding to the options
#'                   \code{"Choudhary P"}, \code{"Escaramis et al."}, \code{"Choudhary NP"} and \code{"Perez-Jaume and Carrasco"} of \code{method}, respectively.
#' @param loa.col colour to be used in the drawing of the LoA. If \code{loa} is set to \code{FALSE}, this argument is ignored. \cr
#'                The default value is \code{"#c27d38"}.
#' @param loess.col colour to be used in the drawing of the loess smooth curve. If \code{loess} is set to \code{FALSE}, this argument is ignored. \cr
#'                The default value is \code{"#cd2c35"}.
#' @param loess.span smoothness parameter for \code{\link[stats]{loess.smooth}}. \cr
#'                   The default value is 2/3.
#' @param legend logical indicating whether a legend should be added outside the plot. If all \code{tdi}, \code{ub} and \code{loa}
#'               are set to \code{FALSE}, this argument is ignored. \cr
#'               The default value is \code{FALSE}.
#' @param inset specifies how far the legend is inset from the plot margins (to be passed to \code{inset} argument in \code{\link[graphics]{legend}}). \cr
#'              The default value is \code{c(-0.25, 0)}, recommended for 24'' screens with default plot window. For 13'' screens, \code{c(-0.34, 0)} is recommended.
#' @param main overall title for the plot (to be passed to \code{main} argument in \code{\link[base]{plot}}). \cr
#'             The default value is \code{"Bland-Altman plot"}.
#' @param xlab a label for the x-axis (to be passed to \code{xlab} argument in \code{\link[base]{plot}}). \cr
#'             The default value is \code{"Mean"}.
#' @param ylab a label for the y-axis (to be passed to \code{ylab} argument in \code{\link[base]{plot}}). \cr
#'             The default value is \code{"Difference"}.
#' @param xlim the x limits of the plot (to be passed to \code{xlim} argument in \code{\link[base]{plot}}). \cr
#'             The default value, \code{NULL}, indicates that the range of the mean values should be used.
#' @param ylim the y limits of the plot (to be passed to \code{ylim} argument in \code{\link[base]{plot}}). \cr
#'             The default value, \code{NULL}, indicates that the range of the differences values should be used.
#' @param ... other graphical parameters (to be passed to \code{\link[base]{plot}}).
#'
#' @returns A Bland-Altman plot of the data in \code{x} with a solid black line at differences \eqn{=} 0, with differences
#' considered as first level \eqn{-} second level of the variable \code{met} in the call of the function \code{\link[TDIagree]{TDI}}.
#'
#' @section Note:
#' A call to \code{\link[graphics]{par}} is used in this method. Notice that the arguments
#' \code{font.lab} and \code{las} are always set to 2 and 1 respectively. Moreover,
#' if \code{legend} is \code{TRUE}, \code{mar} is set to \code{c(4, 4, 2, 9)}.
#'
#' @importFrom graphics par abline lines
#' @importFrom stats loess.smooth
#' @importFrom plotfunctions isColor
#'
#' @examples
#' # normal data
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
#' plot(tdi)
#'
#' # enhance plot
#' plot(tdi, xlim = c(20, 70), ylim = c(-20, 30), tdi = TRUE, ub = TRUE,
#'      method = c("es", "pe"), ub.pc = "b_cb", loa = TRUE, loa.col = "red",
#'      legend = TRUE)
#'
#'
#' # non-normal data
#'
#' tdi.aml <- TDI(AMLad, mrd, id, met, rep, p = c(0.85, 0.95), boot.type = "cluster",
#'                dec.est = 4, R = 1000)
#' plot(tdi.aml)
#'
#' # enhance plot
#' plot(tdi.aml, method = c("choudhary p", "pe"), tdi = TRUE, ub = TRUE, legend = TRUE,
#'      main = "Bland-Altman plot of the MRD")
#'
#'
#' @references Altman, D. G., & Bland, J. M. (1983). Measurement in medicine: the analysis of method comparison studies. Journal of the Royal Statistical Society Series D: The Statistician, 32(3), 307-317.
#'
#'    Bland, J. M., & Altman, D. (1986). Statistical methods for assessing agreement between two methods of clinical measurement. The Lancet, 327(8476), 307-310.
#'
#'    Perez‐Jaume, S., & Carrasco, J. L. (2015). A non‐parametric approach to estimate the total deviation index for non‐normal data. Statistics in Medicine, 34(25), 3318-3335.
#'
#' @exportS3Method TDIagree::plot


plot.tdi <- function(x, tdi = FALSE, ub = FALSE, loa = FALSE, method = NULL,
                     ub.pc = NULL, p = NULL, loess = FALSE,
                     method.col = NULL, loa.col = "#c27d38",
                     loess.col = "#cd2c35", loess.span = 2/3,
                     legend = FALSE, inset = c(-0.24, 0), main = "Bland-Altman plot",
                     xlab = "Mean", ylab = "Difference", xlim = NULL, ylim = NULL, ...) {

  # argument type checks
  if(!inherits(x, "tdi")) {
    stop("'x' must be of class 'tdi'")
  }
  if (!is.logical(tdi)) {
    stop("'tdi' must be logical")
  }
  if (!is.logical(ub)) {
    stop("'ub' must be logical")
  }
  if (!is.logical(loa)) {
    stop("'loa' must be logical")
  }
  if (is.null(method)) { method <- x$params$method }
  choices.method <- tolower(c("Choudhary P", "Escaramis et al.", "Choudhary NP", "Perez-Jaume and Carrasco"))
  method <- tolower(method)
  method <- match.arg(method, choices.method, several.ok = TRUE)
  if (!is.logical(loess)) {
    stop("'loess' must be logical")
  }
  if (!is.null(method.col)) {
    if (any(!isColor(method.col))) {
      stop("'method.col' must be NULL or only contain valid colours")
    }
  }
  if (length(loa.col) > 1) {
    stop("'loa.col' must be of length 1")
    if (!isColor(loa.col)) {
      stop("'loa.col' must only contain a valid colour")
    }
  }
  if (length(loess.col) > 1) {
    stop("'loess.col' must be of length 1")
    if (!isColor(loess.col)) {
      stop("'loess.col' must only contain a valid colour")
    }
  }
  if (!is.numeric(loess.span)) {
    stop("'loess.span' must be numeric")
  }
  if (!is.null(p) & !is.numeric(p)) {
    stop("'p' must be NULL or numeric")
  }
  if (!is.null(ub.pc) & !is.character(ub.pc)) {
    stop("'ub.pc' must be NULL or character")
  }
  if (!is.character(main)) {
    stop("'main' must be a character")
  }
  if (!is.character(xlab)) {
    stop("'xlab' must be a character")
  }
  if (!is.character(ylab)) {
    stop("'ylab' must be a character")
  }
  if (!is.null(xlim)) {
    if (any(is.na(xlim)) | length(xlim) != 2) {
      stop("'xlim' must be NULL or a 2-dimensional vector containing no NAs")
    }
  }
  if (!is.null(ylim)) {
    if (any(is.na(ylim)) | length(ylim) != 2) {
      stop("'ylim' must be NULL or a 2-dimensional vector containing no NAs")
    }
  }

  # reset par at exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # get variables
  result <- x$result
  alpha <- x$params$alpha

  # warning checks
  if (ub & !x$params$ub){
    warning("'ub' set to TRUE but no UBs estimated in the TDI() call, proceeding assuming ub = FALSE")
    ub <- FALSE
  }

  if (!tdi & !ub & !loa & !loess & legend) {
    warning("legend asked but no estimates added to the plot, proceeding assuming legend = FALSE")
    legend <- FALSE
  }

  # stop checks
  if (tdi | ub) {
    if(is.null(p)) {
      p <- result$p[1]
    } else {
      if (length(p) > 1) {
        stop("only one value of 'p' supported")
      }
      if (!p %in% result$p > 0) {
        p.aux <- paste(result$p, collapse = " ")
        stop(paste0("value of 'p' not used in the TDI() call. Possible value(s) of 'p': ", p.aux))
      }
    }
    if ("choudhary p" %in% method) {
      if (!"tdi.ch.p" %in% names(result)) {
        stop("parametric method of Choudhary not estimated in the TDI() call. Rerun without 'choudhary p' in 'method'")
      }
    }
    if ("escaramis et al." %in% method) {
      if (!"tdi.es.p" %in% names(result)) {
        stop("method of Escaramis et al. not estimated in the TDI() call. Rerun without 'escaramis et al.' in 'method'")
      }
    }
    if ("choudhary np" %in% method) {
      if (!"tdi.ch.np" %in% names(result)) {
        stop("non-parametric method of Choudhary not estimated in the TDI() call. Rerun without 'choudhary np' in 'method'")
      }
    }
    if ("perez-jaume and carrasco" %in% method) {
      if (!"tdi.pc.np" %in% names(result)) {
        stop("method of Perez-Jaume and Carrasco not estimated in the TDI() call. Rerun without 'perez-jaume and carrasco' in 'method'")
      }
      if (ub) {
        ub.pc.options <- NULL
        if ("ub_p_db.pc.np" %in% names(result)) {
          ub.pc.options <- c(ub.pc.options, "p_db", "n_db", "e_db", "b_db")
        }
        if ("ub_p_cb.pc.np" %in% names(result)) {
          ub.pc.options <- c(ub.pc.options, "p_cb", "n_cb", "e_cb", "b_cb")
        }
        if (is.null(ub.pc)) {
          ub.pc <- ub.pc.options[1]
        } else if (!ub.pc %in% ub.pc.options) {
          ub.pc.options.aux <- paste(ub.pc.options, collapse = " ")
          stop(paste0("'ub.pc' must be one of the valid values: ", ub.pc.options.aux))
        }
        ub.pc <- paste0("ub_", ub.pc, ".pc.np")
        if (!ub.pc %in% names(result)) {
          stop("technique asked in 'ub.pc' not estimated in the TDI() call")
        }
      }
    }
    if (is.null(method.col)) {
      colours <- c("#f3df6c", "#9c964a", "#f4b5bd", "#85d4e3")
      names(colours) <- c("choudhary p", "escaramis et al.", "choudhary np", "perez-jaume and carrasco")
    } else{
      if (length(method.col) != length(method)) {
        stop("number of colours provided in 'method.col' not equal to number of methods to be added in the plot")
      }
      colours <- method.col
      names(colours) <- method
    }
  }

  data <- x$data.wide

  nrep <- (ncol(data)-1)/2

  if(nrep == 1){
    d <- data$y_met1 - data$y_met2
    m <- (data$y_met1 + data$y_met2)/2
  } else{
    for (k in 1:nrep) {
      for (l in 1:nrep) {
        dif <- data[, paste0("y_met1rep", k)] - data[, paste0("y_met2rep", l)]
        me <- (data[, paste0("y_met1rep", k)] + data[, paste0("y_met2rep", l)])/2
        if(k == 1 & l == 1) {
          d <- dif
          m <- me
        } else{
          d <- append(d, dif)
          m <- append(m, me)
        }
      }
    }
  }

  if (is.null(xlim)) {
    min.x <- min(m)
    max.x <- max(m)
  }
  else {
    min.x <- xlim[1]
    max.x <- xlim[2]
  }
  if (is.null(ylim)) {
    min.y.v <- c(min(d))
    max.y.v <- c(max(d))
    if (tdi) {
      if ("choudhary p" %in% method) {
        min.y.v <- c(min.y.v, -abs(result[result$p == p, "tdi.ch.p"]))
        max.y.v <- c(max.y.v, abs(result[result$p == p, "tdi.ch.p"]))
      }
      if ("escaramis et al." %in% method) {
        min.y.v <- c(min.y.v, -abs(result[result$p == p, "tdi.es.p"]))
        max.y.v <- c(max.y.v, abs(result[result$p == p, "tdi.es.p"]))
      }
      if ("choudhary np" %in% method) {
        min.y.v <- c(min.y.v, -abs(result[result$p == p, "tdi.ch.np"]))
        max.y.v <- c(max.y.v, abs(result[result$p == p, "tdi.ch.np"]))
      }
      if ("perez-jaume and carrasco" %in% method) {
        min.y.v <- c(min.y.v, -abs(result[result$p == p, "tdi.pc.np"]))
        max.y.v <- c(max.y.v, abs(result[result$p == p, "tdi.pc.np"]))
      }
    }
    if (ub) {
      if ("choudhary p" %in% method) {
        min.y.v <- c(min.y.v, -abs(result[result$p == p, "ub.ch.p"]))
        max.y.v <- c(max.y.v, abs(result[result$p == p, "ub.ch.p"]))
      }
      if ("escaramis et al." %in% method) {
        min.y.v <- c(min.y.v, -abs(result[result$p == p, "ub.es.p"]))
        max.y.v <- c(max.y.v, abs(result[result$p == p, "ub.es.p"]))
      }
      if ("choudhary np" %in% method) {
        min.y.v <- c(min.y.v, -abs(result[result$p == p, "ub.ch.np"]))
        max.y.v <- c(max.y.v, abs(result[result$p == p, "ub.ch.np"]))
      }
      if ("perez-jaume and carrasco" %in% method) {
        min.y.v <- c(min.y.v, -abs(result[result$p == p, ub.pc]))
        max.y.v <- c(max.y.v, abs(result[result$p == p, ub.pc]))
      }
    }
    if (loa) {
      z <- qnorm(1-alpha/2, 0, 1)
      min.y.v <- c(min.y.v, mean(d) + z*sd(d), mean(d) - z*sd(d))
      max.y.v <- c(max.y.v, mean(d) + z*sd(d), mean(d) - z*sd(d))
    }
    min.y <- min(min.y.v)
    max.y <- max(max.y.v)
  }
  else {
    min.y <- ylim[1]
    max.y <- ylim[2]
  }

  if (legend) {
    legend.legend <- legend.col <- legend.lty <- c()

    par(font.lab = 2, las = 1, mar = c(4, 4, 2, 9))
    plot(m, d, main = main, xlab = xlab, ylab = ylab,
         xlim = c(min.x, max.x), ylim = c(min.y, max.y), ...)
    abline(h = 0)

    # add horizontal solid lines for TDIs
    if (tdi) {
      if ("choudhary p" %in% method) {
        abline(h = result[result$p == p, "tdi.ch.p"], col = colours["choudhary p"])
        abline(h = -result[result$p == p, "tdi.ch.p"], col = colours["choudhary p"])
        legend.legend <- c(legend.legend, "Ch.P TDI")
        legend.col <- c(legend.col, colours["choudhary p"])
        legend.lty <- c(legend.lty, 1)
      }
      if ("escaramis et al." %in% method) {
        abline(h = result[result$p == p, "tdi.es.p"], col = colours["escaramis et al."])
        abline(h = -result[result$p == p, "tdi.es.p"], col = colours["escaramis et al."])
        legend.legend <- c(legend.legend, "Es. et al. TDI")
        legend.col <- c(legend.col, colours["escaramis et al."])
        legend.lty <- c(legend.lty, 1)
      }
      if ("choudhary np" %in% method) {
        abline(h = result[result$p == p, "tdi.ch.np"], col = colours["choudhary np"])
        abline(h = -result[result$p == p, "tdi.ch.np"], col = colours["choudhary np"])
        legend.legend <- c(legend.legend, "Ch.NP TDI")
        legend.col <- c(legend.col, colours["choudhary np"])
        legend.lty <- c(legend.lty, 1)
      }
      if ("perez-jaume and carrasco" %in% method) {
        abline(h = result[result$p == p, "tdi.pc.np"], col = colours["perez-jaume and carrasco"])
        abline(h = -result[result$p == p, "tdi.pc.np"], col = colours["perez-jaume and carrasco"])
        legend.legend <- c(legend.legend, "P-J&C TDI")
        legend.col <- c(legend.col, colours["perez-jaume and carrasco"])
        legend.lty <- c(legend.lty, 1)
      }
    }

    # add horizontal dashed lines for UB
    if (ub) {
      if ("choudhary p" %in% method) {
        abline(h = result[result$p == p, "ub.ch.p"], col = colours["choudhary p"], lty = 2)
        abline(h = -result[result$p == p, "ub.ch.p"], col = colours["choudhary p"], lty = 2)
        legend.legend <- c(legend.legend, "Ch.P UB")
        legend.col <- c(legend.col, colours["choudhary p"])
        legend.lty <- c(legend.lty, 2)
      }
      if ("escaramis et al." %in% method) {
        abline(h = result[result$p == p, "ub.es.p"], col = colours["escaramis et al."], lty = 2)
        abline(h = -result[result$p == p, "ub.es.p"], col = colours["escaramis et al."], lty = 2)
        legend.legend <- c(legend.legend, "Es. et al. UB")
        legend.col <- c(legend.col, colours["escaramis et al."])
        legend.lty <- c(legend.lty, 2)
      }
      if ("choudhary np" %in% method) {
        abline(h = result[result$p == p, "ub.ch.np"], col = colours["choudhary np"], lty = 2)
        abline(h = -result[result$p == p, "ub.ch.np"], col = colours["choudhary np"], lty = 2)
        legend.legend <- c(legend.legend, "Ch.NP UB")
        legend.col <- c(legend.col, colours["choudhary np"])
        legend.lty <- c(legend.lty, 2)
      }
      if ("perez-jaume and carrasco" %in% method) {
        abline(h = result[result$p == p, ub.pc], col = colours["perez-jaume and carrasco"], lty = 2)
        abline(h = -result[result$p == p, ub.pc], col = colours["perez-jaume and carrasco"], lty = 2)
        legend.legend <- c(legend.legend, "P-J&C UB")
        legend.col <- c(legend.col, colours["perez-jaume and carrasco"])
        legend.lty <- c(legend.lty, 2)
      }
    }

    # add horizontal dotted lines for LoA
    if (loa) {
      z <- qnorm(1-alpha/2, 0, 1)
      abline(h = mean(d) + z*sd(d), col = loa.col, lty = 3)
      abline(h = mean(d) - z*sd(d), col = loa.col, lty = 3)
      legend.legend <- c(legend.legend, paste0(round((1-alpha)*100, 0), "% LoA"))
      legend.col <- c(legend.col, loa.col)
      legend.lty <- c(legend.lty, 3)
    }

    # add dotdashed curve for loess
    if (loess) {
      lines(loess.smooth(m, d, span = loess.span), col = loess.col, lty = 4)
      legend.legend <- c(legend.legend, "LOESS")
      legend.col <- c(legend.col, loess.col)
      legend.lty <- c(legend.lty, 4)
    }

    legend("topright", legend = legend.legend, col = legend.col, lty = legend.lty,
           xpd = TRUE, inset = inset)
  } else{
    par(font.lab = 2, las = 1)
    plot(m, d, main = main, xlab = xlab, ylab = ylab,
         xlim = c(min.x, max.x), ylim = c(min.y, max.y), ...)

    abline(h = 0)

    # add horizontal solid lines for TDIs
    if (tdi) {
      if ("choudhary p" %in% method) {
        abline(h = result[result$p == p, "tdi.ch.p"], col = colours["choudhary p"])
        abline(h = -result[result$p == p, "tdi.ch.p"], col = colours["choudhary p"])
      }
      if ("escaramis et al." %in% method) {
        abline(h = result[result$p == p, "tdi.es.p"], col = colours["escaramis et al."])
        abline(h = -result[result$p == p, "tdi.es.p"], col = colours["escaramis et al."])
      }
      if ("choudhary np" %in% method) {
        abline(h = result[result$p == p, "tdi.ch.np"], col = colours["choudhary np"])
        abline(h = -result[result$p == p, "tdi.ch.np"], col = colours["choudhary np"])
      }
      if ("perez-jaume and carrasco" %in% method) {
        abline(h = result[result$p == p, "tdi.pc.np"], col = colours["perez-jaume and carrasco"])
        abline(h = -result[result$p == p, "tdi.pc.np"], col = colours["perez-jaume and carrasco"])
      }
    }

    # add horizontal dashed lines for UBs
    if (ub) {
      if ("choudhary p" %in% method) {
        abline(h = result[result$p == p, "ub.ch.p"], col = colours["choudhary p"], lty = 2)
        abline(h = -result[result$p == p, "ub.ch.p"], col = colours["choudhary p"], lty = 2)
      }
      if ("escaramis et al." %in% method) {
        abline(h = result[result$p == p, "ub.es.p"], col = colours["escaramis et al."], lty = 2)
        abline(h = -result[result$p == p, "ub.es.p"], col = colours["escaramis et al."], lty = 2)
      }
      if ("choudhary np" %in% method) {
        abline(h = result[result$p == p, "ub.ch.np"], col = colours["choudhary np"], lty = 2)
        abline(h = -result[result$p == p, "ub.ch.np"], col = colours["choudhary np"], lty = 2)
      }
      if ("perez-jaume and carrasco" %in% method) {
        abline(h = result[result$p == p, ub.pc], col = colours["perez-jaume and carrasco"], lty = 2)
        abline(h = -result[result$p == p, ub.pc], col = colours["perez-jaume and carrasco"], lty = 2)
      }
    }

    # add horizontal dotted lines for LoA
    if (loa) {
      z <- qnorm(1-alpha/2, 0, 1)
      abline(h = mean(d) + z*sd(d), col = loa.col, lty = 3)
      abline(h = mean(d) - z*sd(d), col = loa.col, lty = 3)
    }

    # add horizontal dotdashed curve for loess
    if (loess) {
      lines(loess.smooth(m, d, span = loess.span), col = loess.col, lty = 4)
    }
  }
}
