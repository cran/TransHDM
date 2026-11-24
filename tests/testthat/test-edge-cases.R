analysis_args <- function(data, source = NULL, ...) {
  args <- list(
    target_data = data,
    source_data = source,
    Y = "Y", D = "D", M = "M1", X = "X1",
    transfer = FALSE, verbose = FALSE, use_SIS = FALSE
  )
  modifyList(args, list(...))
}

sis_args <- function(data, source = NULL, ...) {
  args <- analysis_args(data, source, ...)
  args$use_SIS <- NULL
  args
}

test_that("analysis functions reject missing, nonnumeric, and nonfinite data", {
  base <- toy_frame()

  missing <- base
  missing$M1[1] <- NA
  expect_error(
    do.call(TransHDM, analysis_args(missing)),
    "Target data contains missing values in: M1",
    fixed = TRUE
  )

  infinite <- base
  infinite$M1[1] <- Inf
  expect_error(
    do.call(SIS, sis_args(infinite)),
    "Target data contains nonfinite values in: M1",
    fixed = TRUE
  )

  text <- base
  text$M2 <- letters[seq_len(nrow(base))]
  expect_error(
    do.call(TransHDM, analysis_args(text, M = c("M1", "M2"))),
    "Target data contains nonnumeric columns: M2",
    fixed = TRUE
  )
})

test_that("analysis functions reject zero variance, unknown names, and bad tuning", {
  base <- toy_frame()

  flat <- base
  flat$D <- 1
  expect_error(
    do.call(TransHDM, analysis_args(flat)),
    "Target data contains zero-variance predictors: D",
    fixed = TRUE
  )

  flat_y <- base
  flat_y$Y <- 3
  expect_error(
    do.call(SIS, sis_args(flat_y)),
    "outcome 'Y' has zero variance",
    fixed = TRUE
  )

  expect_error(
    do.call(TransHDM, analysis_args(base, M = "NOPE")),
    "Some M variables not found in data: NOPE",
    fixed = TRUE
  )
  expect_error(
    do.call(TransHDM, analysis_args(base, p_cutoff = 2)),
    "'p_cutoff' must be a single number in (0, 1].",
    fixed = TRUE
  )
  expect_error(
    do.call(TransHDM, analysis_args(base, ncore = 1.5)),
    "'ncore' must be a positive integer.",
    fixed = TRUE
  )
  expect_error(
    do.call(SIS, sis_args(base, topN = 0)),
    "'topN' must be NULL or a positive integer.",
    fixed = TRUE
  )
  expect_error(
    do.call(TransHDM, analysis_args(base, transfer = TRUE, source = NULL)),
    "requires at least one source dataset",
    fixed = TRUE
  )
  expect_error(
    source_detection(base[1, ], list(base), Y = "Y", D = "D", M = "M1", X = "X1"),
    "at least 2 rows",
    fixed = TRUE
  )
})

test_that("incompatible source columns are rejected before estimation", {
  base <- toy_frame()
  src <- base
  src$M1 <- NULL
  expect_error(
    do.call(TransHDM, analysis_args(base, source = list(base, src), transfer = TRUE)),
    "Source dataset 2 does not contain the same analysis variables as the target data. Missing: M1",
    fixed = TRUE
  )

  extra <- base
  extra$Z <- 1
  expect_error(
    do.call(SIS, sis_args(base, source = list(base, extra), transfer = TRUE)),
    "same columns",
    fixed = TRUE
  )

  reordered <- base[, rev(colnames(base))]
  expect_message(
    combined <- TransHDM:::.combine_sources(list(base, reordered), verbose = TRUE),
    "Combining 2 source datasets into one source (n = 12).",
    fixed = TRUE
  )
  expect_equal(colnames(combined), colnames(base))
  expect_equal(nrow(combined), 12)
  expect_null(TransHDM:::.combine_sources(NULL))
  expect_identical(TransHDM:::.combine_sources(base), base)
})

test_that("one mediator uses the local FDR formula", {
  infer <- fake_inference("M1")
  input_pvalues <- matrix(
    c(infer$P_alpha_SIS[["M1"]], infer$P_beta_SIS[["M1"]]),
    nrow = 1
  )
  input_pvalues <- input_pvalues + matrix(seq_len(2) / 3 * 1e-10, nrow = 1)
  nullprop <- TransHDM:::null_estimation(input_pvalues)
  pmax_i <- max(input_pvalues)
  fdrcut <- pmax_i * nullprop$alpha01 +
    pmax_i * nullprop$alpha10 +
    pmax_i * pmax_i * nullprop$alpha00

  fit <- joint_test(infer, p_cutoff = 1, verbose = FALSE)

  expect_s3_class(fit, "TransHDM")
  expect_equal(fit$all_mediators$mediator, "M1")
  expect_equal(nrow(fit$all_mediators), 1)
  expect_length(fit$IDE_est, 1)
  expect_true(is.finite(fdrcut))
  if (fdrcut <= 1) {
    expect_equal(fit$contributions$mediator, "M1")
    expect_equal(fit$contributions$fdr_est, fdrcut)
    expect_equal(
      fit$contributions$ab_pv,
      max(infer$P_alpha_SIS[["M1"]], infer$P_beta_SIS[["M1"]])
    )
    expect_true(fit$all_mediators$significant)
  } else {
    expect_equal(nrow(fit$contributions), 0)
    expect_false(fit$all_mediators$significant)
  }
})

test_that("joint testing keeps the FDR selection rule and adds fdr_est", {
  one <- joint_test(fake_inference("M1"), p_cutoff = 1, verbose = FALSE)
  expect_s3_class(one, "TransHDM")
  expect_true("fdr_est" %in% names(one$contributions))
  expect_false("fdr_est" %in% names(one$all_mediators))
  expect_equal(nrow(one$all_mediators), 1)
  if (nrow(one$contributions) > 0) {
    expect_equal(
      one$contributions$ab_pv,
      pmax(one$contributions$alpha_pv, one$contributions$beta_pv)
    )
    expect_true(all(one$contributions$fdr_est <= 1))
  }
  est <- setNames(one$effects$estimate, one$effects$effect)
  expect_equal(est[["total"]], est[["indirect"]] + est[["direct"]])

  two <- joint_test(fake_inference(c("M1", "M2")), p_cutoff = 0.05, verbose = FALSE)
  expect_true("fdr_est" %in% names(two$contributions))
  expect_false(any(c("ab_pv", "fdr_est") %in% names(two$all_mediators)))
  if (nrow(two$contributions) > 0) {
    expect_true(all(two$contributions$fdr_est <= 0.05))
  }

  empty <- fake_inference("M1")
  empty$M_SIS <- character(0)
  expect_error(
    joint_test(empty, verbose = FALSE),
    "No mediator passes screening, so joint testing cannot be performed.",
    fixed = TRUE
  )
  expect_error(
    joint_test(list(), verbose = FALSE),
    "mediation_inference",
    fixed = TRUE
  )
  expect_error(
    joint_test(fake_inference("M1"), p_cutoff = 0, verbose = FALSE),
    "'p_cutoff' must be a single number in (0, 1].",
    fixed = TRUE
  )
})

test_that("null estimation keeps a one-row matrix", {
  np <- TransHDM:::null_estimation(matrix(c(0.2, 0.3), nrow = 1))
  expect_true(all(c("alpha00", "alpha01", "alpha10", "alpha1", "alpha2") %in% names(np)))
  expect_error(TransHDM:::null_estimation(c(0.1, 0.2)), "matrix or data frame", fixed = TRUE)
  expect_error(
    TransHDM:::null_estimation(matrix(0.1, nrow = 2, ncol = 1)),
    "2 column",
    fixed = TRUE
  )
})

test_that("a parallel cluster is stopped without closing other connections", {
  con <- file()
  on.exit(close(con), add = TRUE)
  backend <- foreach::getDoParName()
  workers <- foreach::getDoParWorkers()
  cl <- TransHDM:::.start_cluster(1L)
  expect_true(isOpen(con))
  expect_equal(foreach::getDoParName(), backend)
  expect_equal(foreach::getDoParWorkers(), workers)
  TransHDM:::.stop_cluster(cl)
  expect_true(isOpen(con))
  expect_equal(foreach::getDoParName(), backend)
  expect_equal(foreach::getDoParWorkers(), workers)
  classes <- showConnections()[, "class"]
  expect_false(any(grepl("sockconn", classes)))
})

test_that("screening and inference objects reject an empty mediator set", {
  expect_error(
    mediation_inference(list()),
    "class 'SIS'",
    fixed = TRUE
  )
  expect_error(
    mediation_inference(),
    "Either 'screen_result' or 'target_data' must be provided.",
    fixed = TRUE
  )
  expect_error(
    mediation_inference(target_data = toy_frame(), Y = "Y"),
    "'Y', 'D', 'M', and 'X' must be provided",
    fixed = TRUE
  )

  screen <- fake_sis()
  screen$target_SIS$M1 <- NULL
  expect_error(
    mediation_inference(screen, transfer = FALSE, verbose = FALSE),
    "No mediator passes screening",
    fixed = TRUE
  )
  expect_error(
    mediation_inference(fake_sis(), transfer = TRUE, verbose = FALSE),
    "screening result contains none",
    fixed = TRUE
  )
})
