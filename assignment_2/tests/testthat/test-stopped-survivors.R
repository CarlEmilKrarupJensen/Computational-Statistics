test_that("survivor indexing exactly reproduces stopped simulations", {
  for (theta in c(-0.6, -0.3, -0.05, 0, 0.3)) {
    set.seed(6101)
    original <- estimate_ruin(2000, theta, "stopped", diagnostics = TRUE)
    set.seed(6101)
    improved <- estimate_ruin_survivors(2000, theta, diagnostics = TRUE)
    expect_identical(improved, original)
    expect_identical(improved$steps, as.integer(pmin(improved$ruin_times, 100)))
  }
})

test_that("survivor loop handles every path ruining and no path ruining", {
  # These support choices make both cases deterministic despite random draws.
  for (bounds in list(c(-2, -1), c(1, 2))) {
    set.seed(6102)
    original <- estimate_ruin(20, -0.3, "stopped", initial_capital = 0.5,
                              horizon = 4, lower = bounds[1], upper = bounds[2],
                              diagnostics = TRUE)
    set.seed(6102)
    improved <- estimate_ruin_survivors(20, -0.3, initial_capital = 0.5,
                                       horizon = 4, lower = bounds[1], upper = bounds[2],
                                       diagnostics = TRUE)
    expect_identical(improved, original)
    expect_true(all(improved$steps == if (bounds[1] < 0) 1 else 4))
  }
})

test_that("diagnostics are optional and invalid arguments are rejected", {
  result <- estimate_ruin_survivors(10, -0.3)
  expect_named(result, c("summary", "method", "theta"))
  expect_error(estimate_ruin_survivors(1))
  expect_error(estimate_ruin_survivors(10, theta = Inf))
  expect_error(estimate_ruin_survivors(10, horizon = 0))
})
