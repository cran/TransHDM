#' Print method for SIS objects
#'
#' @param x An object of class \code{"SIS"}.
#' @param ... Further arguments (currently not used).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print SIS
#' @export
print.SIS <- function(x, ...) {
  n_sel <- length(x$M_ID_name_SIS)
  n_tot <- length(x$M)
  has_src <- is.data.frame(x$source_SIS)
  cat(sprintf("SIS: %d/%d mediators retained%s\n",
              n_sel, n_tot, if (has_src) " (transfer)" else ""))
  invisible(x)
}


#' @method summary SIS
#' @export
summary.SIS <- function(object, ...) {

  if (!inherits(object, "SIS")) {
    stop("Object must be of class 'SIS'")
  }

  target_SIS <- object$target_SIS
  source_SIS <- object$source_SIS
  M_sel <- object$M_ID_name_SIS

  res <- list(
    selected_mediators = M_sel,
    n_selected = length(M_sel),
    target_dim = if (is.data.frame(target_SIS))
      c(n = nrow(target_SIS), p = ncol(target_SIS)) else NULL,
    source_dim = if (is.data.frame(source_SIS))
      c(n = nrow(source_SIS), p = ncol(source_SIS)) else NULL,
    transfer = is.data.frame(source_SIS)
  )

  class(res) <- c("summary.SIS", "SIS")
  res
}

#' @method print summary.SIS
#' @export
print.summary.SIS <- function(x, ...) {

  cat("\nSure Independence Screening (SIS) summary\n")
  cat("----------------------------------------\n")

  ## selected mediators
  cat("Selected mediators:\n")
  cat("  Number selected:", x$n_selected %||% 0, "\n")

  if (!is.null(x$selected_mediators) && length(x$selected_mediators) > 0) {
    names_txt <- paste(x$selected_mediators, collapse = ", ")
    cat(
      paste(
        strwrap(
          paste0("  Names: ", names_txt),
          width = getOption("width", 80)
        ),
        collapse = "\n"
      ),
      "\n"
    )
  } else {
    cat("  Names: (none)\n")
  }

  ## target info
  if (!is.null(x$target_dim)) {
    cat("\nTarget data after SIS:\n")
    cat("  Sample size (n):", x$target_dim["n"], "\n")
  }

  ## source info
  if (!is.null(x$source_dim)) {
    cat("\nSource data after SIS:\n")
    cat("  Sample size (n):", x$source_dim["n"], "\n")
    cat("\nTransfer learning: YES\n")
  } else {
    cat("\nTransfer learning: NO\n")
  }

  cat("----------------------------------------\n")

  invisible(x)
}


