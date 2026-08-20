#' Summary of Mediation Inference
#'
#' @param object An object of class \code{"mediation_inference"} returned by
#'   \code{\link{mediation_inference}}.
#' @param top Integer, maximum number of mediators to display (default: 10).
#' @param digits Number of digits for rounding estimates (default: 4).
#' @param ... Further arguments (unused).
#'
#' @return
#' An object of class \code{"summary.mediation_inference"}.
#'
#' @method summary mediation_inference
#' @export
summary.mediation_inference <- function(object, top = 10, digits = 4, ...) {

  if (!inherits(object, "mediation_inference")) {
    stop("'object' must be of class 'mediation_inference'.")
  }

  stopifnot(
    is.list(object),
    all(c("DE_est", "alpha_SIS_est", "beta_SIS_est",
          "P_alpha_SIS", "P_beta_SIS", "M", "M_SIS") %in% names(object))
  )

  # Build a quick alpha*beta table for the screened mediators
  alpha_vec <- object$alpha_SIS_est[object$M_SIS]
  beta_vec  <- object$beta_SIS_est[object$M_SIS]
  p_alpha   <- object$P_alpha_SIS[object$M_SIS]
  p_beta    <- object$P_beta_SIS[object$M_SIS]
  ab_vec    <- alpha_vec * beta_vec

  tab <- data.frame(
    mediator   = object$M_SIS,
    alpha      = alpha_vec,
    alpha_pv   = p_alpha,
    beta       = beta_vec,
    beta_pv    = p_beta,
    alpha_beta = ab_vec,
    row.names  = NULL,
    stringsAsFactors = FALSE
  )

  # Sort by absolute mediation effect
  tab <- tab[order(abs(tab$alpha_beta), decreasing = TRUE), , drop = FALSE]

  n_total <- length(object$M)
  n_screened <- length(object$M_SIS)
  show_n <- min(top, n_screened)

  res <- list(
    DE_est        = object$DE_est,
    n_total       = n_total,
    n_screened    = n_screened,
    top_table     = tab[seq_len(show_n), , drop = FALSE],
    top           = top,
    digits        = digits
  )

  class(res) <- "summary.mediation_inference"
  res
}


#' @method print summary.mediation_inference
#' @export
print.summary.mediation_inference <- function(x, ...) {

  p_cols <- c("alpha_pv", "beta_pv")
  digits <- x$digits

  cat("====================================================\n")
  cat("    Summary of Mediation Estimation and Testing       \n")
  cat("====================================================\n\n")

  cat("Mediators:\n")
  cat("  Total specified: ", x$n_total, "\n")
  cat("  Retained after screening: ", x$n_screened, "\n")

  cat("\nDirect Effect (DE):\n")
  cat("  DE_est =", round(x$DE_est, digits), "\n")

  if (x$n_screened == 0) {
    cat("\n(No mediators retained after screening)\n")
    cat("\n====================================================\n")
    return(invisible(x))
  }

  cat("\nTop", min(x$top, x$n_screened),
      "mediators by |alpha * beta|:\n\n")

  tab_show <- x$top_table

  for (nm in names(tab_show)) {
    if (is.numeric(tab_show[[nm]])) {
      if (nm %in% p_cols) {
        tab_show[[nm]] <- formatC(
          tab_show[[nm]],
          format = "e",
          digits = digits
        )
      } else {
        tab_show[[nm]] <- round(tab_show[[nm]], digits)
      }
    }
  }

  print(tab_show, row.names = FALSE)

  cat("\nNote:\n")
  cat("  alpha   : exposure-mediator effect\n")
  cat("  beta    : mediator-outcome effect\n")
  cat("  P-values are shown in scientific notation\n")

  cat("\n====================================================\n")

  invisible(x)
}


#' Print method for mediation_inference objects
#'
#' @param x An object of class \code{"mediation_inference"}.
#' @param ... Further arguments (currently not used).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print mediation_inference
#' @export
print.mediation_inference <- function(x, ...) {
  n_scr <- length(x$M_SIS)
  # n_tot <- length(x$M)
  de <- x$DE_est
  cat(sprintf("Mediation inference: %d mediators tested | DE = %.4f\n",
              n_scr, de))
  invisible(x)
}
