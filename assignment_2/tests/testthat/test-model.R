test_that("original and tilted densities normalize, including the uniform limit", {
  expect_equal(integrate(uniform_density, -1.9, 2)$value, 1)
  for (theta in c(-2, -0.282, 0, 0.4, 2)) {
    expect_equal(integrate(function(x) tilted_density(x, theta), -1.9, 2)$value, 1)
  }
  x <- seq(-2, 2.1, length.out = 100)
  expect_equal(tilted_density(x, 0), uniform_density(x))
  expect_equal(log_mgf(0), 0)
  expect_equal(tilted_mean(0), 0.05)
  for (theta in c(-0.3, 0, 0.3)) {
    numerical <- (log_mgf(theta + 1e-5) - log_mgf(theta - 1e-5)) / 2e-5
    expect_equal(tilted_mean(theta), numerical, tolerance = 1e-7)
  }
})

test_that("inverse transform agrees with the theoretical mean", {
  set.seed(5101)
  for (theta in c(-0.6, -0.282, 0, 0.4)) {
    x <- r_tilted(100000, theta)
    expect_true(all(x >= -1.9 & x <= 2))
    expect_lt(abs(mean(x) - tilted_mean(theta)), 5 * sd(x) / sqrt(length(x)))
  }
})
