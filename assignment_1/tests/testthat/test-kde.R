test_that("KDE is nonnegative", {
  x <- c(-1, 0, 1)

  fit <- kernel_density_estimator(
    x,
    h = 0.5,
    n = 1000
  )

  expect_true(all(fit$y >= 0))
})

test_that("KDE approximately integrates to one", {
  set.seed(1)

  x <- rnorm(100)

  fit <- kernel_density_estimator(
    x,
    h = 0.5,
    n = 5000
  )

  area <- sum(
    diff(fit$x) *
      (head(fit$y, -1) + tail(fit$y, -1)) / 2
  )

  expect_equal(
    area,
    1,
    tolerance = 1e-3
  )
})

test_that("kde_at agrees with grid estimator", {
  x <- c(-1, 0, 1)
  h <- 0.5

  fit <- kernel_density_estimator(
    x,
    h,
    n = 100
  )

  y_direct <- kde_at(
    eval_points = fit$x,
    x = x,
    h = h
  )

  expect_equal(
    fit$y,
    y_direct
  )
})