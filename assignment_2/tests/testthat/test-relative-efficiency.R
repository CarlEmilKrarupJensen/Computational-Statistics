test_that("relative efficiency balances variance and time against a fixed reference", {
  example <- data.frame(theta = c(-0.4, -0.3), contribution_variance = c(3, 2),
                        cost_per_sample = c(1, 3))
  result <- relative_survivor_efficiency(example, reference_theta = -0.3)
  expect_equal(result$relative_variance, c(1.5, 1))
  expect_equal(result$relative_runtime, c(1/3, 1))
  expect_equal(result$relative_cost, c(0.5, 1))
  expect_equal(result$predicted_time_saved_pct, c(50, 0))
  expect_equal(result$relative_cost,
               result$contribution_variance * result$cost_per_sample / 6)
  expect_error(relative_survivor_efficiency(example, -0.2), "Exactly one")
})
