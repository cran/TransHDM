## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(TransHDM)

## ----data---------------------------------------------------------------------
seed <- 1
set.seed(seed)

# ---------------- Simulation Parameters ---------------- #
p_m <- 50     # num of mediators
n <- 100      # num of target samples
rho <- 0.1    # rho for simulation data generation
p_x <- 5      # num of covariates
n_s <- 300    # num of source samples

# ---------------- Target Data Generation ---------------- #
target_sim <- gen_simData_homo(n = n, p_x = p_x, p_m = p_m, rho = rho)
target_data <- target_sim$data

# true effect of target data
true_effect <- target_sim$coef$beta2 * target_sim$coef$alpha1

# column names
M_col <- paste0("M", 1:p_m)
X_col <- paste0("X", 1:p_x)

# ---------------- Source Data Generation ---------------- #
# source, transferable, homogeneous
s_data <- gen_simData_homo(n = n_s, p_x = p_x, p_m = p_m, rho = rho,
                          source = TRUE, transferable = TRUE, h = 2, seed=seed)$data

# source, not transferable, homogeneous
s_f_data <- gen_simData_homo(n = n_s, p_x = p_x, p_m = p_m, rho = rho,
                            source = TRUE, transferable = FALSE, h = 2, seed=seed)$data

# source, transferable, heterogeneous
s_h_data <- gen_simData_hetero(n = n_s, p_x = p_x, p_m = p_m, rho = rho,
                              source = TRUE, transferable = TRUE, h = 2, seed=seed)$data

# source, not transferable, heterogeneous
s_hf_data <- gen_simData_hetero(n = n_s, p_x = p_x, p_m = p_m, rho = rho,
                               source = TRUE, transferable = FALSE, h = 2, seed=seed)$data

## ----show_data----------------------------------------------------------------
# ---------------- Show Data ---------------- #
# show true mediator effect
true_effect

# # show target data
# head(target_data)

## ----detection----------------------------------------------------------------
detect_all <- source_detection(
  target_data = target_data,
  source_data = list(s_data, s_f_data,s_h_data, s_hf_data),
  Y = "Y",
  D = "D",
  M = M_col,
  X = X_col,
  kfold = 5,
  C0 = 0.01,
  verbose = TRUE
)
summary(detect_all)

## ----TransHDM_notrans---------------------------------------------------------
set.seed(seed)

# mediation analysis without transfer learning
res_n <- TransHDM(
  target_data = target_data,
  source_data = NULL,
  Y = "Y",
  D = "D",
  M = M_col,
  X = X_col,
  transfer = FALSE,
  topN = NULL,
  dblasso_SIS = FALSE,
  verbose = TRUE,
  ncore = 1
)
summary(res_n)

## ----TransHDM_transfer--------------------------------------------------------
# mediation analysis with transfer learning (using homogeneous data)
res_t <- TransHDM(
  target_data = target_data,
  source_data = s_data,
  Y = "Y",
  D = "D",
  M = M_col,
  X = X_col,
  transfer = TRUE,
  topN = NULL,
  dblasso_SIS = FALSE,
  verbose = TRUE,
  ncore = 1
)
summary(res_t)

## ----TransHDM_hetero, eval = FALSE--------------------------------------------
# # mediation analysis with transfer learning (using heterogeneous data)
# res_h <- TransHDM(
#   target_data = target_data,
#   source_data = s_h_data,
#   Y = "Y",
#   D = "D",
#   M = M_col,
#   X = X_col,
#   transfer = TRUE,
#   topN = NULL,
#   dblasso_SIS = FALSE,
#   verbose = TRUE,
#   ncore = 1
# )
# summary(res_h)

## ----TransHDM_paral, eval = FALSE---------------------------------------------
# res_p <- TransHDM(
#   target_data = target_data,
#   source_data = s_data,
#   Y = "Y",
#   D = "D",
#   M = M_col,
#   X = X_col,
#   transfer = TRUE,
#   topN = NULL,
#   dblasso_SIS = FALSE,
#   verbose = TRUE,
#   ncore = 4
# )
# summary(res_h)

## ----SIS, eval = FALSE--------------------------------------------------------
# # SIS without transfer learning
# SIS_n <- SIS(
#   target_data = target_data,
#   source_data = NULL,
#   Y = "Y",
#   D = "D",
#   M = M_col,
#   X = X_col,
#   topN = 10,
#   transfer = FALSE,
#   verbose = TRUE,
#   ncore = 1,
#   dblasso_method = FALSE
# )
# summary(SIS_n)
# 
# # SIS with transfer learning
# SIS_t <- SIS(
#   target_data = target_data,
#   source_data = s_data,
#   Y = "Y",
#   D = "D",
#   M = M_col,
#   X = X_col,
#   topN = 10,
#   transfer = TRUE,
#   verbose = TRUE,
#   ncore = 1,
#   dblasso_method = FALSE
# )
# summary(SIS_t)

## ----data_reg-----------------------------------------------------------------
library(MASS)
n_target <- 1000
n_source <- 2000
p <- 20

Sigma <- 0.2^abs(outer(1:p, 1:p, "-"))  # Autocorrelation structure, weak correlation
X_target <- mvrnorm(n_target, mu = rep(0, p), Sigma = Sigma)
X_source <- mvrnorm(n_source, mu = rep(0, p), Sigma = Sigma)

# Construct signal coefficients
# First 3 variables are strong signals, next 2 are weak signals, rest are zero
beta <- c(1.5, -1, 1.0, 0.5, -0.5, rep(0, p-5))

# Construct response variables with noise
y_target <- X_target %*% beta + rnorm(n_target, sd = 1)
y_source <- X_source %*% beta + rnorm(n_source, sd = 1)

# Build target/source lists
target <- list(x = X_target, y = y_target)
source <- list(x = X_source, y = y_source)

## ----lasso--------------------------------------------------------------------
# Fit lasso without transfer learning
coef_n_l <- lasso(target = target, transfer = FALSE, lambda = 'lambda.1se')
summary(coef_n_l)

# Fit lasso with transfer learning
coef_t_l <- lasso(target = target, source = source, transfer = TRUE, lambda = 'lambda.1se')
summary(coef_t_l)

## ----dblasso------------------------------------------------------------------
# Fit dblasso without transfer learning
coef_n_d <- dblasso(target = target, transfer = FALSE, lambda = 'lambda.1se')
summary(coef_n_d)

# Fit dblasso with transfer learning
coef_t_d <- dblasso(target = target, source = source, transfer = TRUE, lambda = 'lambda.1se')
summary(coef_t_d)

