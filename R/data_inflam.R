#' Inflammatory Mediation Benchmark Data: Target Cohort
#'
#' @description
#' Simulated **target** cohort for a high-dimensional mediation analysis in
#' which an inflammatory biomarker (exposure) influences a disease severity
#' index (outcome) through a panel of candidate molecular mediators, adjusting
#' for clinical covariates.
#'
#' @format A data frame with 50 rows (patients) and 25 columns:
#' \itemize{
#'   \item \code{DSI}: outcome — disease severity index.
#'   \item \code{IB}: exposure — inflammatory biomarker.
#'   \item \code{CYT_01}–\code{SIG_04}: 20 candidate mediators across five
#'     functional groups (cytokines \code{CYT}, immune markers \code{IMM},
#'     oxidative stress \code{OXI}, metabolic \code{MET}, signalling \code{SIG}).
#'   \item \code{BIS}, \code{PRS}, \code{BRS}: 3 covariates (baseline
#'     inflammatory status, physiological reserve score, baseline risk score).
#' }
#'
#' @details
#' **Mediation model:**
#' \deqn{IB \rightarrow Molecular\ Mediators \rightarrow DSI}
#'
#' This is the small-sample target cohort (*n* = 50) that benefits from
#' transfer learning using the larger external cohorts. The mediator columns
#' have been randomly permuted, so the identity of the truly active mediators
#' is hidden; the ground truth is stored in \code{\link{inflam_effect}}.
#'
#' @note
#' **Demonstration data only.** Simulated for methodological demonstration of
#' the \code{TransHDM} workflow, not for drawing biological conclusions.
#'
#' @source
#' Generated with \code{\link{gen_simData_homo}} (n = 50, p_x = 3, p_m = 20,
#' rho = 0.1, seed = 123). See \code{gen_data.R} in the package source.
#'
#' @seealso
#' \code{\link{inflam_external1}}, \code{\link{inflam_external2}} for the
#' external cohorts, and \code{\link{inflam_effect}} for the ground truth.
#'
#' @usage data(inflam_target)
#'
#' @examples
#' data(inflam_target)
#' str(inflam_target)
"inflam_target"


#' Inflammatory Mediation Benchmark Data: Transferable External Cohort
#'
#' @description
#' Simulated **transferable** external cohort for high-dimensional mediation
#' analysis. It shares the mediator–outcome structure with the target cohort,
#' making it suitable for transfer learning.
#'
#' @format A data frame with 200 rows (patients) and 25 columns, with the same
#'   structure as \code{\link{inflam_target}}:
#'   \itemize{
#'     \item \code{DSI}: outcome — disease severity index.
#'     \item \code{IB}: exposure — inflammatory biomarker.
#'     \item \code{CYT_01}–\code{SIG_04}: 20 candidate mediators.
#'     \item \code{BIS}, \code{PRS}, \code{BRS}: 3 covariates.
#'   }
#'
#' @details
#' This large external cohort (*n* = 200) is generated under a transferable
#' setting: its exposure–mediator (\eqn{\alpha}) and mediator–outcome
#' (\eqn{\beta}) coefficients match those of the target for the active
#' mediators. It is intended as the source in \code{\link{TransHDM}} with
#' \code{transfer = TRUE}.
#'
#' @note
#' **Demonstration data only.** Simulated for methodological demonstration of
#' the \code{TransHDM} workflow.
#'
#' @source
#' Generated with \code{\link{gen_simData_homo}} (n = 200, p_x = 3, p_m = 20,
#' rho = 0.1, seed = 123, source = TRUE, transferable = TRUE).
#'
#' @seealso
#' \code{\link{inflam_target}} for the target cohort, \code{\link{inflam_external2}}
#' for the non-transferable external cohort.
#'
#' @usage data(inflam_external1)
#'
#' @examples
#' data(inflam_external1)
#' str(inflam_external1)
"inflam_external1"


#' Inflammatory Mediation Benchmark Data: Non-Transferable External Cohort
#'
#' @description
#' Simulated **non-transferable** external cohort for high-dimensional mediation
#' analysis. Its mediator structure differs from the target cohort, and the
#' covariate distribution is shifted, so it should not be used for transfer.
#'
#' @format A data frame with 100 rows (patients) and 25 columns, with the same
#'   structure as \code{\link{inflam_target}}:
#'   \itemize{
#'     \item \code{DSI}: outcome — disease severity index.
#'     \item \code{IB}: exposure — inflammatory biomarker.
#'     \item \code{CYT_01}–\code{SIG_04}: 20 candidate mediators.
#'     \item \code{BIS}, \code{PRS}, \code{BRS}: 3 covariates.
#'   }
#'
#' @details
#' This external cohort (*n* = 100) is generated under a non-transferable
#' setting with covariate shift. Its exposure–mediator and mediator–outcome
#' coefficients differ from the target, so naive transfer would degrade
#' performance. It is used with \code{\link{source_detection}} to demonstrate
#' detection of non-transferable sources.
#'
#' @note
#' **Demonstration data only.** Simulated for methodological demonstration of
#' the \code{TransHDM} workflow.
#'
#' @source
#' Generated with \code{\link{gen_simData_hetero}} (n = 100, p_x = 3, p_m = 20,
#' rho = 0.1, seed = 999, source = TRUE, transferable = FALSE).
#'
#' @seealso
#' \code{\link{inflam_target}} for the target cohort, \code{\link{inflam_external1}}
#' for the transferable external cohort.
#'
#' @usage data(inflam_external2)
#'
#' @examples
#' data(inflam_external2)
#' str(inflam_external2)
"inflam_external2"


#' Inflammatory Mediation Benchmark Data: Ground Truth Effects
#'
#' @description
#' Ground-truth mediation effects for the target cohort of the inflammatory
#' mediation benchmark dataset. Maps each candidate mediator to its true
#' exposure–mediator effect (\eqn{\alpha}), mediator–outcome effect
#' (\eqn{\beta}), and indirect effect (\eqn{\alpha \times \beta}).
#'
#' @format A data frame with 20 rows (one per mediator) and 4 columns:
#' \itemize{
#'   \item \code{mediator}: mediator name.
#'   \item \code{alpha}: true exposure–mediator effect.
#'   \item \code{beta}: true mediator–outcome effect.
#'   \item \code{alpha_beta}: true indirect (mediation) effect.
#' }
#'
#' @details
#' Four mediators have non-zero \eqn{\alpha \times \beta} and are the true
#' active mediators. Because the mediator columns were randomly permuted during
#' data generation, these active mediators are scattered across the candidate
#' names — as in a real study, where the ground truth is unknown at analysis
#' time. This table lets users benchmark \code{\link{TransHDM}} against the
#' truth.
#'
#' @seealso
#' \code{\link{inflam_target}}, \code{\link{inflam_external1}},
#' \code{\link{inflam_external2}} for the associated cohorts.
#'
#' @usage data(inflam_effect)
#'
#' @examples
#' data(inflam_effect)
#' inflam_effect[inflam_effect$alpha_beta != 0, ]
"inflam_effect"
