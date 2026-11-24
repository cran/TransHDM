test_that("print and summary methods describe each fitted class", {
  fit <- fake_transhdm()
  expect_output(print(fit), "2 mediator")
  expect_output(print(fit), "TE=")
  sm <- summary(fit, top = 1)
  expect_s3_class(sm, "summary.TransHDM")
  expect_equal(nrow(sm$contributions), 1)
  expect_output(print(sm), "fdr_est")
  expect_output(print(sm), "ab_pv")

  empty <- fake_transhdm(0)
  expect_output(print(empty), "0 mediator")
  expect_output(print(summary(empty)), "No mediators selected")

  sis <- fake_sis()
  expect_output(print(sis), "SIS: 1/2 mediators retained")
  expect_output(print(summary(sis)), "Number selected: 1")
  expect_output(print(summary(sis)), "M1")
  expect_output(print(summary(sis)), "Transfer learning: NO")
  expect_output(print(summary(fake_sis(TRUE))), "Transfer learning: YES")
  expect_error(TransHDM:::summary.SIS(list()), "class 'SIS'", fixed = TRUE)

  null_count <- summary(sis)
  null_count$n_selected <- NULL
  null_count$selected_mediators <- NULL
  expect_output(print(null_count), "Number selected: 0")

  infer <- fake_inference("M1")
  expect_output(print(infer), "1 mediators tested")
  expect_output(print(summary(infer)), "Direct Effect")
  expect_error(
    TransHDM:::summary.mediation_inference(list()),
    "mediation_inference",
    fixed = TRUE
  )
})

test_that("plot methods return ggplot objects for every TransHDM display", {
  fit <- fake_transhdm()
  for (type in c("overall", "mediator", "pvalue", "alpha_beta")) {
    p <- plot(fit, type = type, label = "none")
    expect_s3_class(p, "ggplot")
  }
  expect_s3_class(plot(fit, type = "alpha_beta", label = "significant"), "ggplot")
  expect_s3_class(plot(fit, type = "alpha_beta", label = "all"), "ggplot")

  empty <- fake_transhdm(0)
  expect_s3_class(plot(empty, type = "mediator"), "ggplot")
  expect_s3_class(plot(empty, type = "pvalue"), "ggplot")

  broken <- fit
  broken$effects <- broken$effects[1:2, ]
  expect_error(plot(broken, type = "overall"), "effects", fixed = TRUE)
  broken <- fit
  broken$contributions <- NULL
  expect_error(plot(broken, type = "mediator"), "contributions", fixed = TRUE)
  expect_error(plot(broken, type = "pvalue"), "contributions", fixed = TRUE)
  broken <- fit
  broken$all_mediators <- NULL
  expect_error(plot(broken, type = "alpha_beta"), "all_mediators", fixed = TRUE)
  expect_error(plot(fit, type = "not_a_plot"), "should be one of")
})

test_that("lasso and dblasso print and summary report the fitted coefficients", {
  set.seed(5)
  x <- matrix(rnorm(40 * 3), 40, 3, dimnames = list(NULL, c("D", "M1", "X1")))
  y <- x[, "D"] * 0.8 + rnorm(40)
  target <- list(x = x, y = y)
  source <- list(x = x[1:20, , drop = FALSE], y = y[1:20] + 0.2)

  fit_1se <- suppressWarnings(TransHDM:::lasso(target, lambda = "lambda.1se"))
  fit_min <- suppressWarnings(TransHDM:::lasso(target, lambda = "lambda.min"))
  fit_transfer <- suppressWarnings(TransHDM:::lasso(
    target, source, transfer = TRUE, lambda = "lambda.1se"
  ))

  expect_s3_class(fit_1se, "lasso")
  expect_s3_class(fit_min, "lasso")
  expect_s3_class(fit_transfer, "lasso")
  expect_named(fit_1se, c("(Intercept)", "D", "M1", "X1"))
  expect_true(all(is.finite(fit_1se)))
  expect_true(all(is.finite(fit_transfer)))
  expect_output(print(fit_1se), "LASSO:")
  expect_output(print(summary(fit_1se)), "Lasso regression summary")
  expect_error(TransHDM:::summary.lasso(c(1, 2)), "Invalid lasso object")

  db_1se <- suppressWarnings(TransHDM:::dblasso(target, lambda = "lambda.1se"))
  db_min <- suppressWarnings(TransHDM:::dblasso(target, lambda = "lambda.min"))
  expect_s3_class(db_1se, "dblasso")
  expect_s3_class(db_min, "dblasso")
  expect_named(db_1se$dbcoef.hat, c("(Intercept)", "D", "M1", "X1"))
  expect_true(all(db_1se$P.value >= 0 & db_1se$P.value <= 1))
  expect_true(all(db_1se$se.est >= 0))
  expect_true(all(db_1se$CI$lb <= db_1se$CI$ub))
  expect_output(print(db_1se), "Debiased LASSO:")
  expect_output(print(summary(db_1se)), "Debiased Lasso Inference Summary")
  expect_error(TransHDM:::summary.dblasso(list()), "class 'dblasso'", fixed = TRUE)
})
