#' Bland-Altman plot
#'
#' @description
#' This function creates a Bland-Altman plot, which is used to evaluate the agreement among the quantitative measures taken by two raters.
#' The plot displays the mean of the measurements from both raters in the x-axis and the differences between the measures taken by the two raters in the y-axis.
#'
#' @param x input object of class \code{tdi} resulting from a call of the function \code{\link[TDIagree]{TDI}}.
#' @param main overall title for the plot (to be passed to \code{main} argument in \code{\link[base]{plot}}).
#' @param xlim the \eqn{x} limits of the plot (to be passed to \code{xlim} argument in \code{\link[base]{plot}}). \cr
#'             The default value, NULL, indicates that the range of the mean values should be used.
#' @param ylim the \eqn{y} limits of the plot (to be passed to \code{ylim} argument in \code{\link[base]{plot}}). \cr
#'             The default value, NULL, indicates that the range of the differences values should be used.
#' @param ... other graphical parameters (to be passed to \code{\link[base]{plot}}).
#'
#' @returns A Bland-Altman plot of the data in \code{x}.
#'
#' @section Note:
#' Future versions of the package will include the TDI and UB estimates in the plot.
#'
#' @importFrom graphics abline
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
#' plot(tdi)
#'
#' # enhance plot
#' par(las = 1, font.lab = 2, font.axis = 2)
#' plot(tdi, xlim = c(10, 70), ylim = c(-60, 80), pch = 16, col = "red")
#'
#' @references Altman, D. G., & Bland, J. M. (1983). Measurement in medicine: the analysis of method comparison studies. Journal of the Royal Statistical Society Series D: The Statistician, 32(3), 307-317.
#'
#' @exportS3Method TDIagree::plot


plot.tdi <- function(x, main = NULL, xlim = NULL, ylim = NULL, ...) {

  if(!inherits(x, "tdi")) {
    stop("x must be of class 'tdi'")
  }

  dec.est <- x$params$dec.est
  dec.p <- x$params$dec.p

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
    min.y <- min(d)
    max.y <- max(d)
  }
  else {
    min.y <- ylim[1]
    max.y <- ylim[2]
  }

  if (is.null(main)){
    plot(m, d, main = paste("Bland-Altman plot"), xlab = "Mean", ylab = "Difference",
         xlim = c(min.x, max.x), ylim = c(min.y, max.y), ...)
    abline(h = 0, lty = 2, col = "darkorchid3")
  } else{
    plot(m, d, main = main, xlab = "Mean", ylab = "Difference",
         xlim = c(min.x, max.x), ylim = c(min.y, max.y), ...)
    abline(h = 0, lty = 2, col = "darkorchid3")
  }
}
