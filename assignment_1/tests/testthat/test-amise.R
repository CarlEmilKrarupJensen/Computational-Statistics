test_that("closed-form and numerical R(f'') agree", {
  x <- c(-1, -0.2, 0.5, 1.3)
  r <- 0.7

  numerical <- Rf2_gaussian_numerical(x, r)
  closed <- Rf2_gaussian_closed(x, r)

  expect_equal(
    closed,
    numerical,
    tolerance = 1e-6
  )
})

test_that("AMISE bandwidth is positive and finite", {
  set.seed(1)

  x <- rnorm(100)

  h <- amise_bandwidth(x)

  expect_true(is.finite(h))
  expect_gt(h, 0)
})