# =========================================================================
# Test: TransHDM full pipeline vs stepwise decomposition
# =========================================================================

library(TransHDM)

cat("=== Test: TransHDM pipeline equivalence ===\n\n")

# ---- Setup: generate test data ----
cat("Generating test data ...\n")
set.seed(123)

target_data <- gen_simData_homo(
  n = 50, p_x = 3, p_m = 20, rho = 0.1, seed = 123
)$data

source_data <- gen_simData_homo(
  n = 100, p_x = 3, p_m = 20, rho = 0.1,
  source = TRUE, transferable = TRUE, seed = 123
)$data

M_vars <- paste0("M", 1:20)
X_vars <- paste0("X", 1:3)


# =========================================================================
# Test A: Full TransHDM pipeline (with SIS + transfer)
# =========================================================================
cat("\n--- Test A: Full TransHDM (with SIS + transfer) ---\n")

set.seed(42)
result_full <- TransHDM(
  target_data = target_data,
  source_data = source_data,
  Y = "Y", D = "D", M = M_vars, X = X_vars,
  transfer = TRUE, ncore = 1, topN = 10, verbose = FALSE
)

stopifnot(inherits(result_full, "TransHDM"))
stopifnot(is.data.frame(result_full$contributions))
stopifnot(is.data.frame(result_full$effects))
stopifnot(nrow(result_full$effects) == 4L)
stopifnot(all(c("indirect", "direct", "total", "pe") %in% result_full$effects$effect))
cat("  -> Full pipeline: OK  (", nrow(result_full$contributions), "mediators identified)\n")


# =========================================================================
# Test B: Stepwise pipeline — SIS + step2 + step3 (with transfer)
# =========================================================================
cat("\n--- Test B: SIS + mediation_inference + joint_test (with transfer) ---\n")

set.seed(42)
sis <- SIS(
  target_data = target_data,
  source_data = source_data,
  Y = "Y", D = "D", M = M_vars, X = X_vars,
  transfer = TRUE, topN = 10, verbose = FALSE
)
stopifnot(inherits(sis, "SIS"))
stopifnot(is.data.frame(sis$target_SIS))
stopifnot(length(sis$M_ID_name_SIS) <= 10L)
cat("  -> SIS: OK  (", length(sis$M_ID_name_SIS), "mediators retained)\n")

infer <- mediation_inference(
  screen_result = sis,
  transfer = TRUE, verbose = FALSE, ncore = 1
)
stopifnot(inherits(infer, "mediation_inference"))
stopifnot(is.numeric(infer$DE_est), length(infer$DE_est) == 1L)
stopifnot(length(infer$alpha_SIS_est) == length(M_vars))
stopifnot(length(infer$beta_SIS_est)  == length(M_vars))
stopifnot(identical(names(infer$alpha_SIS_est), M_vars))
cat("  -> mediation_inference: OK\n")

result_stepwise <- joint_test(
  inference_result = infer,
  p_cutoff = 0.05, verbose = FALSE
)

summary(result_stepwise)

stopifnot(inherits(result_stepwise, "TransHDM"))
cat("  -> joint_test: OK  (", nrow(result_stepwise$contributions), "mediators identified)\n")


# =========================================================================
# Test C: Stepwise pipeline — without SIS (use all mediators)
# =========================================================================
cat("\n--- Test C: mediation_inference + joint_test (without SIS) ---\n")

infer_ns <- mediation_inference(
  target_data = target_data,
  source_data = source_data,
  Y = "Y", D = "D", M = M_vars, X = X_vars,
  transfer = TRUE, verbose = FALSE, ncore = 1
)
stopifnot(inherits(infer_ns, "mediation_inference"))
stopifnot(length(infer_ns$M_SIS) == length(M_vars))  # all mediators retained
cat("  -> mediation_inference (no SIS): OK  (", length(infer_ns$M_SIS), "mediators)\n")

result_stepwise_ns <- joint_test(
  inference_result = infer_ns,
  p_cutoff = 0.05, verbose = FALSE
)
stopifnot(inherits(result_stepwise_ns, "TransHDM"))
cat("  -> joint_test (no SIS): OK  (", nrow(result_stepwise_ns$contributions),
    "mediators identified)\n")


# =========================================================================
# Test D: Without transfer (target only)
# =========================================================================
cat("\n--- Test D: Full pipeline without transfer ---\n")

result_notrans <- TransHDM(
  target_data = target_data,
  source_data = NULL,
  Y = "Y", D = "D", M = M_vars, X = X_vars,
  transfer = FALSE, ncore = 1, use_SIS = FALSE, verbose = FALSE
)
stopifnot(inherits(result_notrans, "TransHDM"))
cat("  -> Full pipeline (no transfer, no SIS): OK\n")

# Stepwise without transfer and without SIS
infer_notrans <- mediation_inference(
  target_data = target_data,
  Y = "Y", D = "D", M = M_vars, X = X_vars,
  transfer = FALSE, verbose = FALSE, ncore = 1
)
stopifnot(inherits(infer_notrans, "mediation_inference"))
result_notrans_sw <- joint_test(infer_notrans, p_cutoff = 0.05, verbose = FALSE)
stopifnot(inherits(result_notrans_sw, "TransHDM"))
cat("  -> Stepwise (no transfer, no SIS): OK\n")


# =========================================================================
# Test E: Equivalence — full vs stepwise (with SIS + transfer)
# =========================================================================
cat("\n--- Test E: Equivalence check (full vs stepwise) ---\n")

# effects table
stopifnot(isTRUE(all.equal(result_full$effects, result_stepwise$effects,
                           tolerance = 1e-10)))
cat("  -> effects: MATCH\n")

# contributions
stopifnot(isTRUE(all.equal(result_full$contributions, result_stepwise$contributions,
                           tolerance = 1e-10)))
cat("  -> contributions: MATCH\n")

# IDE estimates
stopifnot(isTRUE(all.equal(result_full$IDE_est, result_stepwise$IDE_est,
                           tolerance = 1e-10)))
cat("  -> IDE_est: MATCH\n")

# DE estimate
stopifnot(isTRUE(all.equal(result_full$DE_est, result_stepwise$DE_est,
                           tolerance = 1e-10)))
cat("  -> DE_est: MATCH\n")


# =========================================================================
# Summary
# =========================================================================
cat("\n========================================\n")
cat("All tests passed.\n")
cat("Summary:\n")
cat("  Full TransHDM pipeline:              OK\n")
cat("  SIS -> mediation_inference -> joint_test (transfer):    OK\n")
cat("  mediation_inference -> joint_test (no SIS):             OK\n")
cat("  mediation_inference -> joint_test (no transfer/SIS):    OK\n")
cat("  Equivalence full == stepwise:        OK\n")
cat("========================================\n")
