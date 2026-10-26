library(testthat)

test_that("the C++ gradient agrees with the R gradient", {
  # test_dir runs this file from tests/project_B, so go up two folders.
  Rcpp::sourceCpp("../../project_B/gradient.cpp")
  problem <- simulate_doses(100, "continuous", seed = 12)$problem
  theta <- c(0.4, 0.8, 0.1, 1.4)
  index <- c(2L, 5L, 9L)

  # Both functions must use the SAME parameters and observation indices.
  r_gradient <- gradient(theta, problem, index)
  cpp_gradient <- dose_gradient_cpp(theta, problem$z, problem$y, index)
  expect_equal(cpp_gradient, r_gradient, tolerance = 1e-10)
})
