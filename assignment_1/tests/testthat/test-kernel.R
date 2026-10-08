library(testthat)

test_that("Epanechnikov kernel has correct support", {
  expect_equal(
    epanechnikov_kernel(c(-2, -1.1, 1.1, 2)),
    c(0, 0, 0, 0)
  )
})

test_that("Epanechnikov kernel is symmetric", {
  u <- seq(-1, 1, length.out = 101)

  expect_equal(
    epanechnikov_kernel(u),
    epanechnikov_kernel(-u)
  )
})

test_that("Epanechnikov kernel integrates to one", {
  val <- integrate(
    epanechnikov_kernel,
    lower = -1,
    upper = 1
  )$value

  expect_equal(val, 1, tolerance = 1e-8)
})

test_that("Epanechnikov constants are correct", {
  expect_equal(RK_epanechnikov(), 3 / 5)
  expect_equal(mu2_epanechnikov(), 1 / 5)
})