#' Sure Independence Screening for High-Dimensional Mediation Analysis
#'
#' @description
#' **Step 1** of the \code{\link{TransHDM}} pipeline: dimension reduction via
#' Sure Independence Screening (SIS). Mediators are ranked based on the product
#' of their marginal associations with the exposure (alpha) and the outcome (beta).
#' Both associations are lasso point estimates. The top-ranked mediators are
#' retained for downstream de-biased lasso estimation
#' (\code{\link{mediation_inference}}, Step 2) and joint testing
#' (\code{\link{joint_test}}, Step 3).
#'
#' The function supports transfer learning, allowing information from one or more
#' source datasets to be leveraged to improve screening stability and robustness
#' in the target dataset. A list of source datasets is row-bound into one source
#' before screening.
#'
#' @param target_data A data frame containing the target dataset. All variables must be numeric.
#' @param source_data A data frame, or a list of data frames, containing source
#'   datasets (optional, default: NULL). All variables must be numeric and use
#'   the same column names as \code{target_data}. A list is row-bound into one
#'   source dataset before screening.
#' @param Y Character string specifying the outcome variable name.
#' @param D Character string specifying the exposure (treatment) variable name.
#' @param M Character vector specifying mediator variable names.
#' @param X Character vector specifying covariate variable names.
#' @param topN An integer specifying the number of mediators to retain after screening.
#' If \code{NULL}, the number is automatically determined as
#' \eqn{\lceil 2n / \log(n) \rceil}, where \eqn{n} is the target sample size.
#' @param transfer A logical value (default: FALSE) indicating whether to apply transfer
#' learning by incorporating the source dataset in the screening procedure.
#' @param verbose A logical value (default: TRUE) controlling whether progress messages
#' are printed to the console.
#' @param ncore An integer (default: 1) specifying the number of CPU cores for parallel
#' computation.
#'
#' @return A list of class \code{"SIS"} with the following components:
#' \itemize{
#'   \item \code{target_SIS}: A data frame containing the outcome, exposure, selected
#'   mediators, and covariates from the target dataset.
#'   \item \code{source_SIS}: A data frame containing the same variables from the
#'   pooled source dataset if \code{transfer = TRUE}; otherwise \code{NULL}.
#'   \item \code{M_ID_name_SIS}: A character vector of selected mediator names.
#'   \item \code{Y}: The resolved outcome variable name.
#'   \item \code{D}: The resolved exposure variable name.
#'   \item \code{M}: The full mediator variable names (before screening).
#'   \item \code{X}: The resolved covariate variable names.
#' }
#'
#' @seealso \code{\link{TransHDM}}, \code{\link{mediation_inference}},
#'   \code{\link{joint_test}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' # Target data
#' M_target <- matrix(rnorm(200 * 50), nrow = 200)
#' colnames(M_target) <- paste0("M", 1:50)
#'
#' target_data <- data.frame(
#'   Y = rnorm(200),
#'   D = rnorm(200),
#'   M_target,
#'   X1 = rnorm(200)
#' )
#'
#' # Source data
#' M_source <- matrix(rnorm(300 * 50), nrow = 300)
#' colnames(M_source) <- paste0("M", 1:50)
#'
#' source_data <- data.frame(
#'   Y = rnorm(300),
#'   D = rnorm(300),
#'   M_source,
#'   X1 = rnorm(300)
#' )
#'
#' # Run SIS
#' result <- SIS(
#'   target_data = target_data,
#'   source_data = source_data,
#'   Y = "Y",
#'   D = "D",
#'   M = paste0("M", 1:50),
#'   X = "X1",
#'   transfer = TRUE,
#'   topN = 10
#' )
#'
#' result$M_ID_name_SIS
#' }
#'
#' @export


SIS <- function(
    target_data,
    source_data = NULL,
    Y,
    D,
    M,
    X,
    topN = NULL,
    transfer = FALSE,
    verbose = TRUE,
    ncore = 1
){

  checked <- .check_analysis_input(
    target_data = target_data,
    source_data = source_data,
    Y = Y, D = D, M = M, X = X,
    transfer = transfer,
    ncore = ncore,
    topN = topN,
    verbose = verbose
  )
  Y <- checked$Y
  D <- checked$D
  M <- checked$M
  X <- checked$X

  source_data <- .combine_sources(
    source_data,
    verbose = verbose && isTRUE(transfer)
  )

  if (verbose)
    message("Step 1: Sure Independence Screening ...  (",
            format(Sys.time(), "%X"), ")")

  if (ncore > 1) {
    cl <- .start_cluster(ncore)
    on.exit(.stop_cluster(cl), add = TRUE)
  }


  ## ---------- resolve columns ----------
  .resolve_cols <- function(data, cols, name) {
    if (is.null(cols)) return(NULL)
    if (is.numeric(cols)) return(colnames(data)[cols])
    if (is.character(cols)) return(cols)
    stop(paste(name, "must be names or indices"))
  }

  Y <- .resolve_cols(target_data, Y, "Y")
  D <- .resolve_cols(target_data, D, "D")
  M <- .resolve_cols(target_data, M, "M")
  X <- .resolve_cols(target_data, X, "X")

  stopifnot(length(Y) == 1, length(D) == 1, length(M) >= 1)

  p_m <- length(M)
  p_x <- ifelse(is.null(X), 0, length(X))
  n_t <- nrow(target_data)

  ## ---------- determine d_0 ----------
  if (is.null(topN)) {
    d_0 <- ceiling(2 * n_t / log(n_t))
  } else {
    d_0 <- topN
  }
  d_0 <- min(p_m, d_0)
  if (d_0 < 1) {
    stop("No mediator passes screening. Increase 'topN' or check the mediator list.")
  }

  ## ---------- extract matrices ----------
  Y_t <- target_data[[Y]]
  D_t <- target_data[[D]]
  M_t <- as.matrix(target_data[, M, drop = FALSE])
  X_t <- if (p_x > 0) as.matrix(target_data[, X, drop = FALSE]) else NULL
  # X_t <- if (p_x > 0) .make_numeric_matrix(target_data, X) else NULL


  if (transfer) {
    Y_s <- source_data[[Y]]
    D_s <- source_data[[D]]
    M_s <- as.matrix(source_data[, M, drop = FALSE])
    X_s <- if (p_x > 0) as.matrix(source_data[, X, drop = FALSE]) else NULL
    # X_s <- if (p_x > 0) .make_numeric_matrix(source_data, X) else NULL
  }

  ## ============================================================
  ## 1. beta_SIS: mediator -> outcome
  ## ============================================================
  beta_SIS <- numeric(p_m)
  names(beta_SIS) <- M

  if (ncore == 1) {

    for (i in seq_len(p_m)) {

      x_t <- cbind(M_t[, i], D_t, X_t)
      colnames(x_t)[1] <- M[i]

      if (!transfer) {
        fit <- lasso(target = list(x = x_t, y = Y_t))
      } else {
        x_s <- cbind(M_s[, i], D_s, X_s)
        colnames(x_s)[1] <- M[i]
        fit <- lasso(
          target = list(x = x_t, y = Y_t),
          source = list(x = x_s, y = Y_s),
          transfer = TRUE
        )
      }
      beta_SIS[i] <- fit[M[i]]
    }

  } else {

    beta_SIS <- unlist(parallel::parLapply(
      cl, seq_len(p_m), .sis_beta_i,
      M_t = M_t, D_t = D_t, X_t = X_t, M = M, Y_t = Y_t,
      transfer = transfer,
      M_s = if (transfer) M_s else NULL,
      D_s = if (transfer) D_s else NULL,
      X_s = if (transfer) X_s else NULL,
      Y_s = if (transfer) Y_s else NULL
    ), use.names = FALSE)
    names(beta_SIS) <- M
  }


  ## ============================================================
  ## 2. alpha_SIS: exposure -> mediator
  ## ============================================================
  alpha_SIS <- numeric(p_m)
  names(alpha_SIS) <- M

  DX_t <- cbind(D_t, X_t)
  colnames(DX_t)[1] <- D

  if (transfer) DX_s <- cbind(D_s, X_s)

  if (ncore == 1) {

    for (i in seq_len(p_m)) {

      if (!transfer) {
        fit <- lasso(target = list(x = DX_t, y = M_t[, i]))
      } else {
        fit <- lasso(
          target = list(x = DX_t, y = M_t[, i]),
          source = list(x = DX_s, y = M_s[, i]),
          transfer = TRUE
        )
      }
      alpha_SIS[i] <- fit[D]
    }

  } else {

    alpha_SIS <- unlist(parallel::parLapply(
      cl, seq_len(p_m), .sis_alpha_i,
      DX_t = DX_t, M_t = M_t, D = D,
      transfer = transfer,
      DX_s = if (transfer) DX_s else NULL,
      M_s = if (transfer) M_s else NULL
    ), use.names = FALSE)
    names(alpha_SIS) <- M
  }


  ## ============================================================
  ## 3. SIS selection
  ## ============================================================
  ab_SIS <- alpha_SIS * beta_SIS
  ID_SIS <- order(abs(ab_SIS), decreasing = TRUE)[seq_len(d_0)]
  M_ID_name_SIS <- M[ID_SIS]
  if (length(M_ID_name_SIS) == 0) {
    stop("No mediator passes screening. Increase 'topN' or check the mediator list.")
  }

  target_SIS <- target_data[, c(Y, D, M_ID_name_SIS, X), drop = FALSE]
  source_SIS <- if (transfer)
    source_data[, c(Y, D, M_ID_name_SIS, X), drop = FALSE]
  else NULL

  if (verbose)
    message("Top ", length(M_ID_name_SIS),
            " mediators selected: ",
            paste(M_ID_name_SIS, collapse = ", "),
            "  (", format(Sys.time(), "%X"), ")")

  res <- list(
    target_SIS = target_SIS,
    source_SIS = source_SIS,
    M_ID_name_SIS = M_ID_name_SIS,
    Y = Y,
    D = D,
    M = M,
    X = X
  )
  class(res) <- "SIS"
  return(res)
}

.sis_beta_i <- function(i, M_t, D_t, X_t, M, Y_t, transfer,
                        M_s = NULL, D_s = NULL, X_s = NULL, Y_s = NULL) {
  x_t <- cbind(M_t[, i], D_t, X_t)
  colnames(x_t)[1] <- M[i]
  if (!transfer) {
    fit <- lasso(target = list(x = x_t, y = Y_t))
  } else {
    x_s <- cbind(M_s[, i], D_s, X_s)
    colnames(x_s)[1] <- M[i]
    fit <- lasso(
      target = list(x = x_t, y = Y_t),
      source = list(x = x_s, y = Y_s),
      transfer = TRUE
    )
  }
  unname(fit[M[i]])
}

.sis_alpha_i <- function(i, DX_t, M_t, D, transfer,
                         DX_s = NULL, M_s = NULL) {
  if (!transfer) {
    fit <- lasso(target = list(x = DX_t, y = M_t[, i]))
  } else {
    fit <- lasso(
      target = list(x = DX_t, y = M_t[, i]),
      source = list(x = DX_s, y = M_s[, i]),
      transfer = TRUE
    )
  }
  unname(fit[D])
}

