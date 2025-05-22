#' Printing tdi objects
#' @description
#' A nice gt table containing the values computed with the function TDI.
#'
#' @param x input object of class \code{tdi} resulting from a call of the function \code{\link[TDIagree]{TDI}}.
#' @param ... currently not in use
#' @returns A nice \strong{gt} table containing the values computed with the function \code{\link[TDIagree]{TDI}}.
#'          The number of decimals for the estimates and the proportions correspond to the arguments \code{dec.est} and \code{dec.p} of the function \code{\link[TDIagree]{TDI}}, respectively.
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
#' tdi
#'
#'
#' # non-normal data
#'
#' tdi.aml <- TDI(AMLad, mrd, id, met, rep, p = c(0.85, 0.95), boot.type = "cluster",
#'                dec.est = 4, R = 1000)
#' tdi.aml
#'
#'
#' @importFrom gt gt fmt_number tab_spanner md contains ends_with tab_footnote cells_column_labels starts_with cols_label cols_align cols_width px tab_style cell_borders cells_body
#' @importFrom katex katex_html
#'
#' @exportS3Method TDIagree::print

print.tdi <- function(x, ...) {

  if(!inherits(x, "tdi")) {
    stop("x must be of class 'tdi'")
  }

  dec.est <- x$params$dec.est
  dec.p <- x$params$dec.p
  # alpha <- x$params$alpha
  # perc <- round((1-alpha)*100, 0)
  # ub.text <- paste0("$\\text{UB}_{", perc, "\\%}(\\hat{\\kappa}_p)$")


  n <- nrow(x$result)
  if(n > 1){
    white_rows <- 2:n
  } else{
    white_rows <- F
  }

  footnote <- "$p$, proportion; $\\hat{\\kappa}_p$, TDI estimate"

  if (x$params$ub) {
    footnote <- paste0(footnote, "; $\\text{UB}$, upper bound")
    if ("perez-jaume and carrasco" %in% x$params$method) {
      footnote <- paste0(footnote, "; $\\text{p}$, strategy based on percentiles; $\\text{n}$, strategy based on the normal distribution; $\\text{e}$, strategy based on the empirical quantiles; $\\text{BC}_a$, strategy based on the bias-corrected and accelerated method")
    }
  }

  result <- x$result |>
    gt() |>
    fmt_number(columns = -"p", decimals = dec.est) |>
    fmt_number(columns = "p", decimals = dec.p) |>
    tab_spanner(label = md("**Differences bootstrap**"),
                columns = contains("_db")) |>
    tab_spanner(label = md("**Cluster bootstrap**"),
                columns = contains("_cb")) |>
    tab_spanner(label = md("**Choudhary P**"),
                columns = ends_with("ch.p"),
                level = 2) |>
    tab_spanner(label = md("**Escaramis *et al.***"),
                columns = ends_with("es.p"),
                level = 2) |>
    tab_spanner(label = md("**Choudhary NP**"),
                columns = ends_with("ch.np"),
                level = 2) |>
    tab_spanner(label = md("**Perez-Jaume and Carrasco**"),
                columns = ends_with("pc.np"),
                level = 2) |>
    tab_spanner(label = md("**Parametric methods**"),
                columns = ends_with(".p")) |>
    tab_spanner(label = md("**Non-parametric methods**"),
                columns = ends_with(".np")) |>
    # cols_label(p ~ md("$p$"),
    #            starts_with("tdi") ~ md("$\\hat{\\kappa}_p$"),
    #            starts_with("ub") ~ md("$\\text{UB}_{95\\%}(\\kappa}_p)$"),
    #            contains("_p") ~ md("$\\text{UB}_{95\\%}^{\\text{p}}(\\kappa_p)$"),
    #            contains("_n") ~ md("$\\text{UB}_{95\\%}^{\\text{n}}(\\kappa_p)$"),
    #            contains("_e") ~ md("$\\text{UB}_{95\\%}^{\\text{e}}(\\kappa_p)$"),
    #            contains("_b") ~ md("$\\text{UB}_{95\\%}^{\\text{BC}_a}(\\kappa_p)$")) |>
    cols_label(p ~ md("$p$"),
               starts_with("tdi") ~ md("$\\hat{\\kappa}_p$"),
               starts_with("ub") ~ md("$\\text{UB}(\\kappa_p)$"),
               contains("_p") ~ md("$\\text{UB}^{\\text{p}}(\\kappa_p)$"),
               contains("_n") ~ md("$\\text{UB}^{\\text{n}}(\\kappa_p)$"),
               contains("_e") ~ md("$\\text{UB}^{\\text{e}}(\\kappa_p)$"),
               contains("_b") ~ md("$\\text{UB}^{\\text{BC}_a}(\\kappa_p)$")) |>
    cols_align("center") |>
    cols_width(everything() ~ px(105)) |>
    tab_style(style = cell_borders(sides = c("right"),
                                   weight = px(2), color = "lightgrey"),
               locations = cells_body(columns = "p")) |>
    tab_style(style = cell_borders(sides = c("top"),
                                   weight = px(2), color = "white"),
              locations = cells_body(rows = white_rows)) |>
    tab_footnote(md(footnote))

  print(result)
}
