#' Printing tdi objects
#' @description
#' A nice gt table containing the values computed with the function TDI.
#'
#' @param x input object of class \code{tdi} resulting from a call of the function \code{\link[TDIagree]{TDI}}.
#' @param ... currently not in use
#' @returns A nice \strong{gt} table containing the values computed with the function \code{\link[TDIagree]{TDI}}.
#'          The number of decimals for the estimates and the probabilities correspond to the arguments \code{dec.est} and \code{dec.p} of the function \code{\link[TDIagree]{TDI}}, respectively.
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
#' tdi
#'
#'
#' @importFrom gt gt fmt_number tab_spanner md ends_with tab_footnote cells_column_labels starts_with cols_label cols_align cols_width px tab_style cell_borders cells_body
#'
#' @exportS3Method TDIagree::print

print.tdi <- function(x, ...) {

  if(!inherits(x, "tdi")) {
    stop("x must be of class 'tdi'")
  }

  dec.est <- x$params$dec.est
  dec.p <- x$params$dec.p

  n <- nrow(x$result)
  if(n > 1){
    white_rows <- 2:n
  } else{
    white_rows <- F
  }

  result <- x$result |>
    gt() |>
    fmt_number(columns = -"p", decimals = dec.est) |>
    fmt_number(columns = "p", decimals = dec.p) |>
    tab_spanner(label = md("**Lin**"),
                columns = ends_with("lin.p")) |>
    tab_spanner(label = md("**Choudhary P**"),
                columns = ends_with("ch.p")) |>
    tab_spanner(label = md("**Escaramis**"),
                columns = ends_with("es.p")) |>
    tab_spanner(label = md("**Choudhary NP**"),
                columns = ends_with("ch.np")) |>
    tab_spanner(label = md("**Perez-Jaume and Carrasco**"),
                columns = ends_with("pc.np")) |>
    tab_spanner(label = md("**Parametric methods**"),
                columns = ends_with(".p")) |>
    tab_spanner(label = md("**Non-parametric methods**"),
                columns = ends_with(".np")) |>
    tab_footnote(footnote = "Total deviation index estimate",
                 locations = cells_column_labels(columns = starts_with("tdi"))) |>
    tab_footnote(footnote = "95% upper bound",
                 locations = cells_column_labels(columns = starts_with("ub."))) |>
    tab_footnote(footnote = "95% upper bound for the bootstrap non-parametric technique based on percentiles",
                 locations = cells_column_labels(columns = starts_with("ub_p"))) |>
    tab_footnote(footnote = "95% upper bound for the bootstrap non-parametric technique based on the normal distribution",
                 locations = cells_column_labels(columns = starts_with("ub_nlog"))) |>
    cols_label("p" ~ md("*p*"),
               starts_with("tdi") ~ md("*{{:kappa:_p}}*"),
               starts_with("ub") ~ "{{UB_95%}}",
               starts_with("ub_p.pc.np") ~ "{{UB[_95%^P]}}",
               starts_with("ub_nlog.pc.np") ~ "{{UB[_95%^N]}}") |>
    cols_align("center") |>
    cols_width(everything() ~ px(80)) |>
    tab_style(style = cell_borders(sides = c("right"),
                                   weight = px(2), color = "lightgrey"),
               locations = cells_body(columns = "p")) |>
    tab_style(style = cell_borders(sides = c("top"),
                                   weight = px(2), color = "white"),
              locations = cells_body(rows = white_rows))

  print(result)
}
