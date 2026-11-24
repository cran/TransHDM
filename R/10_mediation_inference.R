#' Mediation Inference via De-biased Lasso
#'
#' @description
#' **Step 2** of the \code{\link{TransHDM}} pipeline: de-biased lasso estimation
#' following screening (\code{\link{SIS}}, Step 1). It consists of two
#' sub-steps:
#' \itemize{
#'   \item \strong{Step 2.1 — Outcome model (beta):} estimates mediator-outcome
#'     effects by regressing the outcome \code{Y} on all retained mediators,
#'     the exposure \code{D}, and covariates \code{X}.
#'   \item \strong{Step 2.2 — Mediator models (alpha):} estimates
#'     exposure-mediator effects by regressing each mediator on \code{D}
#'     and \code{X}.
#' }
#' Both steps support transfer learning. A list of source datasets is row-bound
#' into one source before estimation. Results are passed to
#' \code{\link{joint_test}} (Step 3) for multiple testing.
#'
#' @details
#' The variance estimator follows Algorithm 3 of Tian and Feng (2023).
#' For Gaussian linear models, it assumes unit error variance.
#' Under the required source-target heterogeneity conditions, pooled
#' source and target design information is used to estimate the
#' correlation structure, while the debiasing correction is based
#' on target-domain data. This assumption is not automatically
#' verified by the function.
#'
#' @param screen_result An object of class \code{"SIS"} returned by
#'   \code{\link{SIS}}. When provided, \code{Y}, \code{D}, \code{M}, and
#'   \code{X} are automatically extracted from it. If \code{NULL},
#'   \code{target_data} (and optionally \code{source_data}) are used
#'   directly without screening, and \code{Y}, \code{D}, \code{M},
#'   \code{X} must be supplied explicitly.
#' @param target_data A data frame containing the target dataset. Required
#'   when \code{screen_result = NULL}. Ignored otherwise.
#' @param source_data A data frame, or a list of data frames, containing source
#'   datasets (optional, default: \code{NULL}). Used only when
#'   \code{transfer = TRUE} and \code{screen_result = NULL}. A list is
#'   row-bound into one source dataset before estimation.
#' @param Y Character string specifying the outcome variable name. Required
#'   only when \code{screen_result = NULL}; otherwise extracted from
#'   \code{screen_result}.
#' @param D Character string specifying the exposure variable name.
#'   Required only when \code{screen_result = NULL}.
#' @param M Character vector specifying the full set of mediator variable
#'   names (before screening). Required only when \code{screen_result = NULL};
#'   otherwise extracted from \code{screen_result}.
#' @param X Character vector specifying covariate variable names.
#'   Required only when \code{screen_result = NULL}.
#' @param transfer A logical value (default: \code{FALSE}) indicating
#'   whether to enable transfer learning.
#' @param verbose A logical value (default: \code{TRUE}) controlling
#'   whether progress messages are printed.
#' @param ncore An integer (default: 1) specifying the number of CPU cores
#'   for parallel computation of mediator models.
#'
#' @return An object of class \code{"mediation_inference"} — a list with:
#' \itemize{
#'   \item \code{DE_est}: Estimated direct effect of D on Y.
#'   \item \code{alpha_SIS_est}: Named numeric vector of alpha estimates
#'     (exposure-mediator effects) for all mediators in \code{M}.
#'     Non-selected mediators are set to zero.
#'   \item \code{P_alpha_SIS}: Named numeric vector of p-values for alpha.
#'   \item \code{beta_SIS_est}: Named numeric vector of beta estimates
#'     (mediator-outcome effects) for all mediators in \code{M}.
#'     Non-selected mediators are set to zero.
#'   \item \code{P_beta_SIS}: Named numeric vector of p-values for beta.
#'   \item \code{M}: The full mediator name vector.
#'   \item \code{M_SIS}: Mediators retained after screening.
#'   \item \code{D}: The exposure variable name.
#' }
#'
#' @seealso \code{\link{TransHDM}}, \code{\link{joint_test}}, \code{\link{SIS}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' target_data <- gen_simData_homo(n = 50, p_x = 3, p_m = 20, rho = 0.1)$data
#' source_data <- gen_simData_homo(n = 100, p_x = 3, p_m = 20, rho = 0.1,
#'   source = TRUE, transferable = TRUE)$data
#' M_vars <- paste0("M", 1:20)
#' X_vars <- paste0("X", 1:3)
#'
#' # With screening (SIS)
#' sis <- SIS(target_data, source_data, Y = "Y", D = "D", M = M_vars,
#'   X = X_vars, transfer = TRUE, topN = 10)
#' infer <- mediation_inference(sis, transfer = TRUE)
#'
#' # Without screening (use all mediators directly)
#' infer_ns <- mediation_inference(
#'   target_data = target_data, source_data = source_data,
#'   Y = "Y", D = "D", M = M_vars, X = X_vars,
#'   transfer = TRUE)
#' }
#'
#' @export
mediation_inference <- function(
    screen_result = NULL,
    target_data = NULL,
    source_data = NULL,
    Y = NULL,
    D = NULL,
    M = NULL,
    X = NULL,
    transfer = FALSE,
    verbose = TRUE,
    ncore = 1
) {

  # ---- Resolve inputs ----
  if (!is.null(screen_result)) {
    # Screening path: extract data and variable names from SIS result
    if (!inherits(screen_result, "SIS")) {
      stop("'screen_result' must be an object of class 'SIS'.")
    }
    target_SIS <- screen_result$target_SIS
    source_SIS <- screen_result$source_SIS

    # Auto-extract variable names from screen_result
    Y <- screen_result$Y
    D <- screen_result$D
    M <- screen_result$M
    X <- screen_result$X

    if (is.null(Y) || is.null(D) || is.null(M) || is.null(X)) {
      stop("'screen_result' does not contain Y, D, M, or X. ",
           "Please update the SIS object or provide them explicitly.")
    }
    if (isTRUE(transfer) && is.null(source_SIS)) {
      stop("'transfer = TRUE' requires source data, but the screening result contains none.")
    }
    M_present <- intersect(M, colnames(target_SIS))
    if (length(M_present) == 0) {
      stop("No mediator passes screening. Increase 'topN' or check that mediator names are present in the data.")
    }
    .check_analysis_input(
      target_data = target_SIS,
      source_data = source_SIS,
      Y = Y, D = D, M = M_present, X = X,
      transfer = isTRUE(transfer) && !is.null(source_SIS),
      ncore = ncore,
      verbose = verbose
    )
  } else {
    # Direct-data path: use raw data without screening
    if (is.null(target_data)) {
      stop("Either 'screen_result' or 'target_data' must be provided.")
    }
    if (is.null(M) || is.null(Y) || is.null(D) || is.null(X)) {
      stop("'Y', 'D', 'M', and 'X' must be provided when 'target_data' is given.")
    }
    checked <- .check_analysis_input(
      target_data = target_data,
      source_data = source_data,
      Y = Y, D = D, M = M, X = X,
      transfer = transfer,
      ncore = ncore,
      verbose = verbose
    )
    Y <- checked$Y
    D <- checked$D
    M <- checked$M
    X <- checked$X
    target_SIS <- target_data
    source_SIS <- .combine_sources(
      source_data,
      verbose = verbose && isTRUE(transfer)
    )
  }

  # ---- Determine mediators retained after screening ----
  M_SIS <- intersect(M, colnames(target_SIS))
  p_m_SIS <- length(M_SIS)

  if (p_m_SIS == 0) {
    stop("No mediator passes screening. Increase 'topN' or check that mediator names are present in the data.")
  }

  if (verbose && p_m_SIS < length(M)) {
    message("After screening, ", p_m_SIS, " / ", length(M),
            " mediators are retained.")
  }

  if (verbose) {
    message("Step 2: De-biased Lasso Estimates ... (",
            format(Sys.time(), "%X"), ")")
  }

  # ---- 2.1 Outcome model (beta) ----
  outcome_x <- setdiff(colnames(target_SIS), Y)

  target_outcome <- list(
    x = as.matrix(target_SIS[, outcome_x, drop = FALSE]),
    y = target_SIS[[Y]]
  )

  if (transfer && !is.null(source_SIS)) {
    source_outcome <- list(
      x = as.matrix(source_SIS[, outcome_x, drop = FALSE]),
      y = source_SIS[[Y]]
    )
    DBLASSO_fit <- dblasso(
      target = target_outcome,
      source = source_outcome,
      transfer = TRUE
    )
  } else {
    DBLASSO_fit <- dblasso(
      target = target_outcome,
      transfer = FALSE
    )
  }

  DE_est <- DBLASSO_fit$dbcoef.hat[D]

  beta_SIS_est <- setNames(rep(0, length(M)), M)
  P_beta_SIS   <- setNames(rep(1, length(M)), M)

  beta_SIS_est[M_SIS] <- DBLASSO_fit$dbcoef.hat[M_SIS]
  P_beta_SIS[M_SIS]   <- DBLASSO_fit$P.value[M_SIS]

  if (verbose) {
    message("Estimation of mediator-outcome effects in the ",
            "outcome model completed.")
  }

  # ---- 2.2 Mediator models (alpha) ----
  alpha_SIS_est <- setNames(rep(0, length(M)), M)
  P_alpha_SIS   <- setNames(rep(1, length(M)), M)

  mediator_x <- c(D, X)
  mediator_x <- intersect(mediator_x, colnames(target_SIS))

  if (ncore == 1) {
    for (m in M_SIS) {
      target_m <- list(
        x = as.matrix(target_SIS[, mediator_x, drop = FALSE]),
        y = target_SIS[[m]]
      )
      if (transfer && !is.null(source_SIS)) {
        source_m <- list(
          x = as.matrix(source_SIS[, mediator_x, drop = FALSE]),
          y = source_SIS[[m]]
        )
        fit <- dblasso(
          target = target_m,
          source = source_m,
          transfer = TRUE
        )
      } else {
        fit <- dblasso(
          target = target_m,
          transfer = FALSE
        )
      }
      alpha_SIS_est[m] <- fit$dbcoef.hat[D]
      P_alpha_SIS[m]   <- fit$P.value[D]
    }
  } else {
    cl <- .start_cluster(ncore)
    on.exit(.stop_cluster(cl), add = TRUE)

    alpha_res <- do.call(rbind, parallel::parLapply(
      cl, M_SIS, .mediator_alpha_i,
      target_SIS = target_SIS,
      mediator_x = mediator_x,
      D = D,
      transfer = transfer,
      source_SIS = source_SIS
    ))

    P_alpha_SIS[M_SIS]   <- alpha_res[, "p"]
    alpha_SIS_est[M_SIS] <- alpha_res[, "a"]
  }

  if (verbose) {
    message("Estimation of exposure-mediator effects in the ",
            "mediator model completed.")
  }

  # ---- Assemble result ----
  result <- list(
    DE_est        = DE_est,
    alpha_SIS_est = alpha_SIS_est,
    P_alpha_SIS   = P_alpha_SIS,
    beta_SIS_est  = beta_SIS_est,
    P_beta_SIS    = P_beta_SIS,
    M             = M,
    M_SIS         = M_SIS,
    D             = D
  )
  class(result) <- "mediation_inference"
  return(result)
}

.mediator_alpha_i <- function(m, target_SIS, mediator_x, D, transfer,
                              source_SIS = NULL) {
  target_m <- list(
    x = as.matrix(target_SIS[, mediator_x, drop = FALSE]),
    y = target_SIS[[m]]
  )
  if (transfer && !is.null(source_SIS)) {
    source_m <- list(
      x = as.matrix(source_SIS[, mediator_x, drop = FALSE]),
      y = source_SIS[[m]]
    )
    fit <- dblasso(
      target = target_m,
      source = source_m,
      transfer = TRUE
    )
  } else {
    fit <- dblasso(target = target_m, transfer = FALSE)
  }
  c(p = unname(fit$P.value[D]), a = unname(fit$dbcoef.hat[D]))
}
