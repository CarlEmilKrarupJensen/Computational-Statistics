test_that("full and stopped weights agree with deterministic density ratios", {
  x <- c(-1, -1.5, 1, 0.2)
  theta <- -0.3
  for (steps in c(2, 4)) {
    expected <- prod(uniform_density(x[seq_len(steps)]) /
                       tilted_density(x[seq_len(steps)], theta))
    expect_equal(exp(log_likelihood_ratio(sum(x[seq_len(steps)]), steps, theta)), expected)
  }
  expect_equal(log_likelihood_ratio(sum(x), 4, 0), 0)
  # Deterministic two-step quadrature verifies unbiasedness and conditioning:
  # midpoint grid, initial capital 1, ruin possible at either step.
  grid <- seq(-1.9, 2, length.out = 601)
  x1 <- head(grid, -1) + diff(grid)[1] / 2
  pairs <- expand.grid(x1 = x1, x2 = x1)
  ruined1 <- 1 + pairs$x1 <= 0
  ruined <- ruined1 | 1 + pairs$x1 + pairs$x2 <= 0
  full <- as.numeric(ruined) * exp(log_likelihood_ratio(pairs$x1 + pairs$x2, 2, theta))
  stopped <- as.numeric(ruined) * exp(log_likelihood_ratio(
    ifelse(ruined1, pairs$x1, pairs$x1 + pairs$x2), ifelse(ruined1, 1, 2), theta))
  mass <- tilted_density(pairs$x1, theta) * tilted_density(pairs$x2, theta) * diff(grid)[1]^2
  exact_grid <- mean(ruined)
  expect_equal(sum(full * mass), exact_grid, tolerance = 1e-6)
  expect_equal(sum(stopped * mass), exact_grid, tolerance = 2e-6)
  expect_lte(sum(stopped^2 * mass), sum(full^2 * mass))
})

test_that("independent unbiased estimators are consistent within MC uncertainty", {
  set.seed(5301)
  results <- lapply(c("crude", "full", "stopped"), function(method) {
    estimate_ruin(80000, -0.28, method, diagnostics = TRUE)
  })
  for (i in 1:2) for (j in (i + 1):3) {
    a <- results[[i]]$summary
    b <- results[[j]]$summary
    expect_lt(abs(a$estimate - b$estimate), 5 * sqrt(a$standard_error^2 + b$standard_error^2))
  }
  stopped <- results[[3]]
  expect_equal(stopped$steps, pmin(stopped$ruin_times, 100))
  expect_equal(stopped$summary$estimator_variance,
               var(stopped$contributions) / length(stopped$contributions))
  expect_equal(results[[1]]$contributions, as.numeric(is.finite(results[[1]]$ruin_times)))
  expect_true(all(results[[2]]$steps == 100))
})
