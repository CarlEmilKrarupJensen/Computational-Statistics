test_that("first crossing includes equality and ignores subsequent recovery", {
  increments <- rbind(c(-1, -1, 2), c(1, -1, 0), c(-3, 2, 2))
  expect_equal(ruin_times(increments, initial_capital = 2), c(2, Inf, 1))
  expect_true(2 + sum(increments[1, ]) > 0)
})

test_that("path simulation includes initial capital", {
  set.seed(5201)
  result <- simulate_paths(10, horizon = 5)
  expect_equal(dim(result$paths), c(10L, 6L))
  expect_equal(result$paths[, 1], rep(30, 10))
})
