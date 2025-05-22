#' TDI.Perez_Jaume_and_Carrasco
#'
#' @param data.wide name of the dataset in its wide format.
#' @param p A value or vector of the probabilities for the TDI.
#' @param ub logical asking whether the UBs should be computed.
#'           Default is TRUE.
#' @param boot.type name of the bootstrap approach(es) to be used. There are two different options when there are replicates:
#'                  to bootstrap the vector of the within-subject differences ("Differences") or to bootstrap at subject level ("Cluster").
#'                  This is, not all the differences coming from the same subject need to be bootstrapped together in the first one but all
#'                  the measurements from the same subjects have to be bootstrapped together in the second one.
#'                  Default is c("Differences", "Cluster").
#' @param alpha confidence level for inference on the TDI.
#'              Default is 0.05.
#' @param type an integer between 1 and 9 selecting one of the nine quantile algorithms (to be passed to quantile).
#'             Recommended 8 for continuous data and 3 for discrete data.
#'             Default is 8.
#' @param R number of bootstrap replicates (to be passed to boot).
#'          Default is 10000.
#'
#' @importFrom stats quantile qnorm sd
#' @importFrom boot boot
#' @importFrom coxed bca
#'
#' @noRd


# auxiliar function for TDI.Perez_Jaume_and_Carrasco

d_boot <- function(dataset, i, p, type){
  return(quantile(dataset[i], p, type = type))
}

c_boot <- function(dataset, i, p, type){
  data.b <- dataset[i, ]

  for (k in 1:((ncol(data.b)-1)/2)) {
    for (l in 1:((ncol(data.b)-1)/2)) {
      dif <- data.b[, paste0("y_met1rep", k)] - data.b[, paste0("y_met2rep", l)]
      if (k == 1 & l == 1) {
        d <- dif
      }
      else{
        d <- append(d, dif)
      }
    }
  }

  return(quantile(abs(d), p, type = type))
}

TDI.Perez_Jaume_and_Carrasco <- function(data.wide, p, ub = TRUE, boot.type = c("Differences", "Cluster"),
                                         alpha = 0.05, type = 8, R = 10000){

  if (ncol(data.wide) > 3) {
    for (k in 1:((ncol(data.wide)-1)/2)) {
      for (l in 1:((ncol(data.wide)-1)/2)) {
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

  # tdi
  tdi.hat <- quantile(abs(d), p, type = type)

  # ub
  ub_p_db.np <- NA
  ub_n_db.np <- NA
  ub_e_db.np <- NA
  ub_b_db.np <- NA

  ub_p_cb.np <- NA
  ub_n_cb.np <- NA
  ub_e_cb.np <- NA
  ub_b_cb.np <- NA

  if (ub) {
    # differences bootsrtap
    if ("differences" %in% boot.type) {
      db <- boot(abs(d), statistic = d_boot, R = R, stype = "i", p = p, type = type)
      tdi.d.boot <- db$t

      ub_p_db.np <- quantile(tdi.d.boot, 1-alpha, type = type)

      if (any(tdi.d.boot == 0)) {
        warning(paste(sum(tdi.d.boot == 0), "zero(s) found in the bootstrap vector, computation of differences UB.nlog with original scale"))
        ub_n_db.np <- tdi.hat+qnorm(1-alpha)*sd(tdi.d.boot)
      }
      else{
        ln.tdi.d.boot <- log(tdi.d.boot)
        ub_n_db.np <- exp(log(tdi.hat)+qnorm(1-alpha)*sd(ln.tdi.d.boot))
      }

      ub_e_db.np <- 2*tdi.hat - quantile(tdi.d.boot, alpha, type = type)
      ub_b_db.np <- bca(tdi.d.boot, conf.level = 1-2*alpha)[2]
    }
    # cluster bootstrap
    if ("cluster" %in% boot.type) {
      cb <- boot(data.wide, statistic = c_boot, R = R, stype = "i", p = p, type = type)

      tdi.c.boot <- cb$t

      ub_p_cb.np <- quantile(tdi.c.boot, 1-alpha, type = type)
      if (any(tdi.c.boot == 0)) {
        warning(paste(sum(tdi.c.boot == 0), "zero(s) found in the bootstrap vector, computation of clustered UB.nlog with original scale"))
        ub_n_cb.np <- tdi.hat+qnorm(1-alpha)*sd(tdi.c.boot)
      } else{
        ln.tdi.c.boot <- log(tdi.c.boot)
        ub_n_cb.np <- exp(log(tdi.hat)+qnorm(1-alpha)*sd(ln.tdi.c.boot))
      }
      ub_e_cb.np <- 2*tdi.hat - quantile(tdi.c.boot, alpha, type = type)
      ub_b_cb.np <- bca(tdi.c.boot, conf.level = 1-2*alpha)[2]
    }
  }
  res <- c(tdi.hat, ub_p_db.np, ub_n_db.np, ub_e_db.np, ub_b_db.np,
                    ub_p_cb.np, ub_n_cb.np, ub_e_cb.np, ub_b_cb.np)
  return(res)
}
