test_that("one mediator and two sources are estimated without a dimension error", {
  target <- gen_simData_hetero(
    n = 36, p_x = 2, rho = 0, p_m = 8, source = FALSE, seed = 11
  )$data
  s1 <- gen_simData_hetero(
    n = 40, p_x = 2, rho = 0, p_m = 8, h = 0,
    source = TRUE, transferable = TRUE, seed = 21
  )$data
  s2 <- gen_simData_hetero(
    n = 32, p_x = 2, rho = 0, p_m = 8, h = 0,
    source = TRUE, transferable = TRUE, seed = 22
  )$data

  messages <- character()
  fit <- withCallingHandlers(
    suppressWarnings(TransHDM(
      target, list(s1, s2),
      Y = "Y", D = "D", M = "M1", X = "X1",
      transfer = TRUE, use_SIS = FALSE, verbose = TRUE, ncore = 1
    )),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl(
    "Combining 2 source datasets into one source (n = 72).",
    messages,
    fixed = TRUE
  )))

  expect_s3_class(fit, "TransHDM")
  expect_true("fdr_est" %in% names(fit$contributions))
  expect_false("fdr_est" %in% names(fit$all_mediators))
  expect_equal(fit$all_mediators$mediator, "M1")
  if (nrow(fit$contributions) > 0) {
    expect_equal(
      fit$contributions$ab_pv,
      pmax(fit$contributions$alpha_pv, fit$contributions$beta_pv)
    )
    expect_true(all(fit$contributions$fdr_est <= 0.05))
  }
  est <- setNames(fit$effects$estimate, fit$effects$effect)
  expect_equal(unname(est[["total"]]), unname(est[["indirect"]] + est[["direct"]]))
  expect_true(is.finite(est[["total"]]))
  expect_equal(nrow(fit$contributions), 0)
  expect_s3_class(plot(fit, type = "overall"), "ggplot")
  expect_s3_class(plot(fit, type = "mediator"), "ggplot")
  expect_output(print(summary(fit)), "No mediators selected")
})

test_that("SIS retains the requested mediators from the supplied list", {
  sim <- gen_simData_hetero(
    n = 40, p_x = 2, rho = 0, p_m = 8, source = FALSE, seed = 8
  )
  mediators <- paste0("M", 1:4)
  res <- suppressWarnings(SIS(
    sim$data,
    Y = "Y", D = "D", M = mediators, X = "X1",
    topN = 2, transfer = FALSE, verbose = FALSE, ncore = 1
  ))

  expect_s3_class(res, "SIS")
  expect_length(res$M_ID_name_SIS, 2)
  expect_true(all(res$M_ID_name_SIS %in% mediators))
  expect_null(res$source_SIS)
  expect_equal(colnames(res$target_SIS), c("Y", "D", res$M_ID_name_SIS, "X1"))
  expect_output(print(res), "2/4 mediators retained")
})

test_that("parallel screening stops only its own cluster", {
  sim <- gen_simData_hetero(
    n = 30, p_x = 2, rho = 0, p_m = 8, source = FALSE, seed = 6
  )
  con <- file()
  on.exit(close(con), add = TRUE)
  backend <- foreach::getDoParName()
  workers <- foreach::getDoParWorkers()
  res <- suppressWarnings(SIS(
    sim$data,
    Y = "Y", D = "D", M = c("M1", "M2"), X = "X1",
    topN = 1, transfer = FALSE, verbose = FALSE, ncore = 2
  ))
  expect_s3_class(res, "SIS")
  expect_length(res$M_ID_name_SIS, 1)
  expect_true(isOpen(con))
  expect_equal(foreach::getDoParName(), backend)
  expect_equal(foreach::getDoParWorkers(), workers)
  classes <- showConnections()[, "class"]
  expect_false(any(grepl("sockconn", classes)))
})
