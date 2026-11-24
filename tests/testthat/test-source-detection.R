.bad_source_data <- function(seed) {
  set.seed(seed)
  n <- 50
  D <- rnorm(n)
  X1 <- rnorm(n)
  M1 <- 0.8 * D + rnorm(n)
  Y <- 1.2 * D + 0.7 * M1 + rnorm(n)
  target <- data.frame(Y, D, M1, X1)

  ns <- 80
  Ds <- rnorm(ns)
  X1s <- rnorm(ns)
  M1s <- rnorm(ns)
  Ys <- -4 * Ds + rnorm(ns, sd = 3)
  source <- data.frame(Y = Ys, D = Ds, M1 = M1s, X1 = X1s)
  list(target = target, source = source)
}

test_that("source_detection follows the loss threshold and reports no transfer", {
  dat <- .bad_source_data(3)
  expect_message(
    res <- suppressWarnings(source_detection(
      dat$target, list(dat$source),
      Y = "Y", D = "D", M = "M1", X = "X1",
      kfold = 2, t = 0, verbose = TRUE
    )),
    "No transferable source is identified.",
    fixed = TRUE
  )

  expect_s3_class(res, "source_detection")
  expect_length(res$source.loss, 1)
  expect_equal(res$T_index, res$source.loss - res$target.valid.loss)
  expect_equal(res$threshold, 0)
  expect_equal(res$transfer.source.id, which(res$T_index <= res$threshold))
  expect_length(res$transfer.source.id, 0)
  expect_equal(dim(res$loss.cv), c(2, 2))
})

test_that("source_detection names a transferable source when the index is below the threshold", {
  dat <- .bad_source_data(1)
  expect_message(
    res <- suppressWarnings(source_detection(
      dat$target, list(dat$source),
      Y = "Y", D = "D", M = "M1", X = "X1",
      kfold = 2, t = 0, verbose = TRUE
    )),
    "Transfer candidate sources: 1",
    fixed = TRUE
  )
  expect_equal(res$transfer.source.id, 1)
  expect_true(res$T_index[1] <= res$threshold)
})

test_that("source_detection accepts one data frame and a list of sources", {
  dat <- .bad_source_data(1)
  extra <- dat$source
  extra[] <- lapply(extra, function(z) z + 0.1)

  one <- suppressWarnings(source_detection(
    dat$target, dat$source,
    Y = "Y", D = "D", M = "M1", X = "X1",
    kfold = 2, t = 0.02, verbose = FALSE
  ))
  two <- suppressWarnings(source_detection(
    dat$target, list(dat$source, extra),
    Y = "Y", D = "D", M = "M1", X = "X1",
    kfold = 2, t = 0.02, verbose = FALSE
  ))

  expect_s3_class(one, "source_detection")
  expect_length(one$source.loss, 1)
  expect_length(two$source.loss, 2)
  expect_equal(two$transfer.source.id, which(two$T_index <= two$threshold))
  expect_equal(two$threshold, 25 / 2 * log1p(0.02))
  expect_equal(dim(two$loss.cv), c(2, 3))
})

test_that("source detection decision is unchanged when every outcome is rescaled", {
  dat <- .bad_source_data(11)
  fit_scaled <- function(scale) {
    target <- dat$target
    source <- dat$source
    target$Y <- target$Y * scale
    source$Y <- source$Y * scale
    set.seed(11)
    suppressWarnings(source_detection(
      target, list(source),
      Y = "Y", D = "D", M = "M1", X = "X1",
      kfold = 2, t = 0.02, verbose = FALSE
    ))
  }

  base <- fit_scaled(1)
  scaled <- fit_scaled(10)

  expect_equal(scaled$threshold, base$threshold)
  expect_equal(scaled$T_index, base$T_index, tolerance = 1e-6)
  expect_equal(scaled$transfer.source.id, base$transfer.source.id)
  expect_equal(base$threshold, 25 / 2 * log1p(0.02))
})

test_that("source_detection rejects missing sources and invalid tuning", {
  dat <- toy_frame()
  expect_error(
    source_detection(dat, NULL, Y = "Y", D = "D", M = "M1", X = "X1"),
    "non-empty list of source data frames",
    fixed = TRUE
  )
  expect_error(
    source_detection(dat, list(), Y = "Y", D = "D", M = "M1", X = "X1"),
    "non-empty list of data frames",
    fixed = TRUE
  )

  src <- dat
  src$M1 <- NULL
  expect_error(
    source_detection(dat, list(src), Y = "Y", D = "D", M = "M1", X = "X1", kfold = 2),
    "Missing: M1",
    fixed = TRUE
  )
  expect_error(
    source_detection(dat, list(dat), Y = "Y", D = "D", M = "M1", X = "X1", kfold = 1),
    "kfold",
    fixed = TRUE
  )
  expect_error(
    source_detection(dat, list(dat), Y = "Y", D = "D", M = "M1", X = "X1", t = -1),
    "'t' must be a single nonnegative number.",
    fixed = TRUE
  )
})

test_that("source detection print, summary, and plot cover empty and nonempty results", {
  found <- fake_source(1L)
  expect_output(print(found), "1/2 sources transferable")
  expect_output(print(summary(found)), "Transferable")
  expect_output(print(summary(found)), "YES")

  none <- found
  none$transfer.source.id <- integer(0)
  expect_output(print(none), "0/2 sources transferable")
  expect_output(print(summary(none)), "Number of transferable sources: 0")
  expect_output(print(summary(none)), "NO")

  empty <- found
  empty$source.loss <- numeric(0)
  empty$T_index <- numeric(0)
  empty$transfer.source.id <- integer(0)
  expect_output(print(summary(empty)), "No source data provided.")

  p_found <- plot(found)
  p_empty <- plot(empty)
  expect_s3_class(p_found, "ggplot")
  expect_s3_class(p_empty, "ggplot")
  expect_error(
    TransHDM:::plot.source_detection(list()),
    "source_detection",
    fixed = TRUE
  )
})
