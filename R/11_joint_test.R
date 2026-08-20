#' Joint Multiple Testing for Mediation Analysis
#'
#' @description
#' **Step 3** (final step) of the \code{\link{TransHDM}} pipeline: joint
#' multiple testing on the alpha and beta p-values obtained from
#' \code{\link{mediation_inference}} (Step 2). It estimates the null
#' proportion, controls the false discovery rate via \code{HDMT::fdr_est},
#' identifies significant mediators, and summarises the direct, indirect,
#' and total effects into a \code{"TransHDM"} object.
#'
#' Together with \code{\link{SIS}} (Step 1) and
#' \code{\link{mediation_inference}} (Step 2), this completes the three-step
#' TransHDM workflow.
#'
#' @param inference_result An object of class \code{"mediation_inference"}
#'   returned by \code{\link{mediation_inference}}.
#' @param p_cutoff A numeric value (default: 0.05) specifying the maximum FDR
#'   threshold for identifying significant mediators.
#' @param verbose A logical value (default: \code{TRUE}) controlling
#'   whether progress messages are printed.
#'
#' @return An object of class \code{"TransHDM"} — a list with:
#' \itemize{
#'   \item \code{contributions}: Data frame of identified mediators with
#'     columns \code{mediator}, \code{alpha}, \code{alpha_pv}, \code{beta},
#'     \code{beta_pv}, \code{alpha_beta}, \code{ab_pv}, \code{pa}.
#'   \item \code{effects}: Data frame summarising indirect, direct, total
#'     effects and proportion mediated (\code{pe}).
#'   \item \code{IDE_est}: Named numeric vector of indirect effect estimates
#'     for all mediators (non-selected set to zero).
#'   \item \code{DE_est}: Estimated direct effect.
#'   \item \code{all_mediators}: Data frame of all SIS-screened mediators with
#'     columns \code{mediator}, \code{alpha}, \code{alpha_pv}, \code{beta},
#'     \code{beta_pv}, \code{significant} (logical).
#' }
#'
#' @seealso \code{\link{TransHDM}}, \code{\link{mediation_inference}}, \code{\link{SIS}}
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
#' # Full pipeline
#' sis <- SIS(target_data, source_data, Y = "Y", D = "D", M = M_vars,
#'   X = X_vars, transfer = TRUE, topN = 10)
#' infer <- mediation_inference(sis, M = M_vars, D = "D", X = X_vars,
#'   Y = "Y", transfer = TRUE)
#' result <- joint_test(infer)
#' plot(result)
#' }
#'
#' @export
joint_test <- function(
    inference_result,
    p_cutoff = 0.05,
    verbose = TRUE
) {

  if (!inherits(inference_result, "mediation_inference")) {
    stop("'inference_result' must be an object of class 'mediation_inference'.")
  }

  # ---- Extract from inference result ----
  DE_est        <- inference_result$DE_est
  alpha_SIS_est <- inference_result$alpha_SIS_est
  P_alpha_SIS   <- inference_result$P_alpha_SIS
  beta_SIS_est  <- inference_result$beta_SIS_est
  P_beta_SIS    <- inference_result$P_beta_SIS
  M             <- inference_result$M
  M_SIS         <- inference_result$M_SIS

  p_m <- length(M)

  # ---- Multiple testing ----
  if (verbose) {
    message("Step 3: Multiple-testing procedure ... (",
            format(Sys.time(), "%X"), ")")
  }

  PA <- cbind(P_alpha_SIS[M_SIS], P_beta_SIS[M_SIS])
  P_value <- apply(PA, 1, max)

  N0 <- dim(PA)[1] * dim(PA)[2]
  input_pvalues <- PA + matrix(runif(N0, 0, 10^{-10}), dim(PA)[1], 2)

  nullprop <- null_estimation(input_pvalues)

  fdrcut <- HDMT::fdr_est(
    nullprop$alpha00,
    nullprop$alpha01,
    nullprop$alpha10,
    nullprop$alpha1,
    nullprop$alpha2,
    input_pvalues,
    exact = 0
  )

  ID_fdr <- which(fdrcut <= p_cutoff)
  selected_M <- M_SIS[ID_fdr]

  alpha_hat_est <- alpha_SIS_est[selected_M]
  beta_hat_est  <- beta_SIS_est[selected_M]
  ab_pv <- P_value[ID_fdr]

  # ---- Effects ----
  alpha_beta_hat_est <- alpha_hat_est * beta_hat_est

  IDE_est <- setNames(rep(0, p_m), M)
  IDE_est[selected_M] <- alpha_beta_hat_est

  TE_est <- sum(IDE_est) + DE_est
  PE_est <- sum(IDE_est) / TE_est

  effects <- data.frame(
    effect   = c("indirect", "direct", "total", "pe"),
    estimate = c(sum(IDE_est), DE_est, TE_est, PE_est),
    stringsAsFactors = FALSE
  )

  contributions <- data.frame(
    mediator   = selected_M,
    alpha      = alpha_hat_est,
    alpha_pv   = P_alpha_SIS[selected_M],
    beta       = beta_hat_est,
    beta_pv    = P_beta_SIS[selected_M],
    alpha_beta = alpha_beta_hat_est,
    ab_pv      = ab_pv,
    pa         = alpha_beta_hat_est / TE_est,
    row.names  = NULL,
    stringsAsFactors = FALSE
  )

  # All SIS-screened mediators (significant + non-significant)
  all_mediators <- data.frame(
    mediator    = M_SIS,
    alpha       = alpha_SIS_est[M_SIS],
    alpha_pv    = P_alpha_SIS[M_SIS],
    beta        = beta_SIS_est[M_SIS],
    beta_pv     = P_beta_SIS[M_SIS],
    significant = M_SIS %in% selected_M,
    row.names   = NULL,
    stringsAsFactors = FALSE
  )

  if (verbose && length(selected_M) > 0) {
    message("Identified mediator(s): ",
            paste(selected_M, collapse = ", "))
  }

  # ---- Return TransHDM object ----
  fit <- list(
    contributions = contributions,
    effects       = effects,
    IDE_est       = IDE_est,
    DE_est        = DE_est,
    all_mediators = all_mediators
  )
  class(fit) <- "TransHDM"
  return(fit)
}

