test_that("gen_simData_hetero returns the documented layout and effect identity", {
  n <- 20
  p_x <- 2
  p_m <- 10
  sim <- gen_simData_hetero(
    n = n, p_x = p_x, rho = 0, p_m = p_m,
    source = FALSE, seed = 7
  )

  expect_s3_class(sim$data, "data.frame")
  expect_equal(dim(sim$data), c(n, 2 + p_m + p_x))
  expect_equal(
    colnames(sim$data),
    c("Y", "D", paste0("M", seq_len(p_m)), paste0("X", seq_len(p_x)))
  )
  expect_true(all(vapply(sim$data, is.numeric, logical(1))))
  expect_true(all(is.finite(as.matrix(sim$data))))

  effect <- sim$effect
  expect_equal(effect$total, effect$direct + effect$indirect)
  expect_equal(effect$indirect, sum(effect$alpha * effect$beta))
  expect_equal(effect$mediation, effect$alpha * effect$beta)
  expect_length(effect$alpha, p_m)
  expect_length(effect$beta, p_m)
  expect_equal(effect$direct, 1)
})

test_that("gen_simData_hetero is reproducible and distinguishes source designs", {
  args <- list(n = 20, p_x = 2, rho = 0.3, p_m = 10, h = 2, seed = 9)
  a <- do.call(gen_simData_hetero, args)
  b <- do.call(gen_simData_hetero, args)
  expect_equal(a, b)

  target <- gen_simData_hetero(
    n = 20, p_x = 2, p_m = 10, source = FALSE, seed = 4
  )
  transferable <- gen_simData_hetero(
    n = 20, p_x = 2, p_m = 10, h = 2,
    source = TRUE, transferable = TRUE, seed = 4
  )
  shifted <- gen_simData_hetero(
    n = 20, p_x = 2, p_m = 10,
    source = TRUE, transferable = FALSE, seed = 4
  )

  expect_equal(transferable$effect$alpha[1:8], target$effect$alpha[1:8])
  expect_equal(transferable$effect$beta[1:8], target$effect$beta[1:8])
  expect_equal(sum(transferable$effect$alpha[9:10] != 0), 2)
  expect_false(identical(shifted$effect$alpha[1:8], target$effect$alpha[1:8]))
  expect_true(all(shifted$effect$alpha[1:8] == 0.5))
  expect_false(identical(transferable$data, target$data))
  expect_true(all(is.finite(as.matrix(transferable$data))))
  expect_true(all(is.finite(as.matrix(shifted$data))))
})

test_that("gen_simData_hetero rejects a mediator dimension below the design", {
  expect_error(
    gen_simData_hetero(n = 20, p_x = 2, p_m = 7, seed = 1)
  )
  expect_error(
    gen_simData_hetero(
      n = 20, p_x = 2, p_m = 10, h = 3,
      source = TRUE, transferable = TRUE, seed = 1
    )
  )
})

test_that("gen_simData_homo matches its effect identity and column layout", {
  sim <- gen_simData_homo(
    n = 16, p_x = 2, rho = 0, p_m = 8, source = FALSE, seed = 3
  )
  expect_equal(dim(sim$data), c(16, 12))
  expect_equal(sim$effect$total, sim$effect$direct + sim$effect$indirect)
  expect_equal(
    gen_simData_homo(n = 16, p_x = 2, p_m = 8, seed = 3),
    gen_simData_homo(n = 16, p_x = 2, p_m = 8, seed = 3)
  )
})
