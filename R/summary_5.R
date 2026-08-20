#' Print method for dblasso objects
#'
#' @param x An object of class \code{"dblasso"}.
#' @param ... Further arguments (currently not used).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print dblasso
#' @export
print.dblasso <- function(x, ...) {
  sig <- sum(x$P.value < 0.05, na.rm = TRUE)
  cat(sprintf("Debiased LASSO: %d significant / %d coefficients (p < 0.05)\n",
              sig, length(x$P.value)))
  invisible(x)
}


#' Summary of Debiased Lasso Inference
#'
#' @param object An object of class \code{"dblasso"}.
#' @param ... Further arguments (currently not used).
#'
#' @return An object of class \code{"summary.dblasso"}.
#'
#' @method summary dblasso
#' @export
summary.dblasso <- function(object, ...) {

  if (!inherits(object, "dblasso")) {
    stop("Object must be of class 'dblasso'")
  }

  est <- object$dbcoef.hat
  se  <- object$se.est
  z   <- est / se
  p   <- object$P.value
  CI  <- object$CI

  tab <- data.frame(
    Estimate  = est,
    Std.Error = se,
    Z.value   = z,
    Pr...z..  = p,
    CI.Lower  = CI$lb,
    CI.Upper  = CI$ub,
    row.names = names(est)
  )

  res <- list(
    coef_table = tab,
    n_coef = length(est)
  )

  class(res) <- "summary.dblasso"
  res
}


#' Print method for summary.dblasso objects
#'
#' @description
#' Displays the coefficient table from a debiased Lasso model summary.
#' Supports fixed decimal places and scientific notation formatting.
#'
#' @param x An object of class \code{"summary.dblasso"}.
#' @param digits Number of digits to display (default is 3).
#' @param fixed Logical. If \code{TRUE}, numbers are printed with fixed decimal
#' places (no scientific notation). If \code{FALSE}, numbers are printed in
#' scientific notation with the specified digits.
#' @param ... Further arguments (currently not used).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print summary.dblasso
#' @export
print.summary.dblasso <- function(x, digits = 3, fixed = FALSE, ...) {
  cat("Debiased Lasso Inference Summary\n")
  cat(strrep("-", 32), "\n", sep = "")
  tab <- x$coef_table
  tab_fmt <- as.data.frame(
    lapply(tab, function(col) {
      if (is.numeric(col)) {
        if (fixed) {
          formatC(col, format = "f", digits = digits)
        } else {
          formatC(col, format = "e", digits = digits)
        }
      } else col
    }),
    row.names = rownames(tab),
    check.names = FALSE
  )
  print(noquote(tab_fmt), quote = FALSE)
  invisible(x)
}
