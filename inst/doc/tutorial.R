## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(TransHDM)

data(inflam_target)
data(inflam_external1)
data(inflam_external2)
data(inflam_effect)

target    <- inflam_target
source    <- inflam_external1
source_nt <- inflam_external2

outcome_var  <- "DSI"    # disease severity index (outcome)
exposure_var <- "IB"     # inflammatory biomarker (exposure)
X_vars       <- c("BIS", "PRS", "BRS")
M_vars       <- setdiff(colnames(target),
                        c(outcome_var, exposure_var, X_vars))

## ground truth effect (target cohort)
inflam_effect

## -----------------------------------------------------------------------------
set.seed(123)
detect <- source_detection(target_data = target,
                           source_data = list(source, source_nt),
                           Y = outcome_var, D = exposure_var, M = M_vars, X = X_vars)
print(detect)
summary(detect)

## ----fig.width=6, fig.height=3.5----------------------------------------------
plot(detect)

## -----------------------------------------------------------------------------
# With transfer learning
set.seed(123)
result_tl <- TransHDM(
  target_data = target,
  source_data = source,
  Y = outcome_var, D = exposure_var,
  M = M_vars, X = X_vars,
  transfer = TRUE,
  topN = 10, p_cutoff = 0.01
)
print(result_tl)
# summary(result_tl)

## -----------------------------------------------------------------------------
# Without transfer learning (target only)
set.seed(123)
result_nt <- TransHDM(
  target_data = target,
  Y = outcome_var, D = exposure_var,
  M = M_vars, X = X_vars,
  transfer = FALSE,
  topN = 10, p_cutoff = 0.01
)
print(result_nt)
summary(result_nt)

## ----fig.width=6, fig.height=3.5----------------------------------------------
# Mediator-wise effect plot
plot(result_tl, type = "mediator")

## ----fig.width=6, fig.height=3.5----------------------------------------------
# P-value visualization
plot(result_tl, type = "pvalue")

## ----fig.width=6, fig.height=3.5----------------------------------------------
# Overall effect decomposition
plot(result_tl, type = "overall")

## ----fig.width=6, fig.height=3.5----------------------------------------------
# Alpha-beta effect plot
plot(result_tl, type = "alpha_beta")

## -----------------------------------------------------------------------------
set.seed(123)
sis <- SIS(
  target_data = target, source_data = source,
  Y = outcome_var, D = exposure_var, M = M_vars, X = X_vars,
  transfer = TRUE, topN = 10, ncore = 1
)
print(sis)
summary(sis)

## -----------------------------------------------------------------------------
med_effect <- mediation_inference(
  screen_result = sis,
  transfer = TRUE
)
print(med_effect)
summary(med_effect, top = 10)

## ----eval=FALSE---------------------------------------------------------------
# mediation_inference(
#   source_data = source, target_data = target,
#   Y = outcome_var, D = exposure_var, M = M_vars, X = X_vars,
#   transfer = TRUE, ncore = 1
# )

## -----------------------------------------------------------------------------
jt <- joint_test(
  inference_result = med_effect,
  p_cutoff = 0.05
)
print(jt)
summary(jt)

