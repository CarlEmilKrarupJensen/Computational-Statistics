test_that("make_folds returns valid fold assignments", {
  set.seed(1)

  x <- rnorm(20)

  folds <- make_folds(x, k = 5)

  expect_length(folds, length(x))
  expect_equal(sort(unique(folds)), 1:5)
})

test_that("likelihood CV score is finite", {
  set.seed(1)

  x <- rnorm(30)
  folds <- make_folds(x, 5)

  score <- cv_loglik_score(
    x = x,
    h = 0.5,
    fold_id = folds
  )

  expect_true(is.finite(score))
})

test_that("LSCV score is finite", {
  set.seed(1)

  x <- rnorm(30)

  score <- lscv_score(
    x = x,
    h = 0.5
  )

  expect_true(is.finite(score))
})


test_that("select_bandwidth_grid maximizes correctly", {
  h_grid <- c(0.1, 0.2, 0.3)

  result <- select_bandwidth_grid(
    h_grid,
    score_function = function(h) -(h - 0.2)^2,
    objective = "max"
  )

  expect_equal(
    result$bandwidth,
    0.2
  )
})

test_that("select_bandwidth_grid minimizes correctly", {
  h_grid <- c(0.1, 0.2, 0.3)

  result <- select_bandwidth_grid(
    h_grid,
    score_function = function(h) (h - 0.2)^2,
    objective = "min"
  )

  expect_equal(
    result$bandwidth,
    0.2
  )
})


test_that("repeated_cv correctly maximizes an objective", {

  x <- 1:20
  h_grid <- c(0.1, 0.2, 0.3, 0.4)

  mock_score <- function(x, h, fold_id, ...) {
    -(h - 0.3)^2
  }

  result <- repeated_cv(
    x = x,
    h_grid = h_grid,
    score_function = mock_score,
    objective = "max",
    k = 5,
    repeats = 3
  )

  expect_equal(
    result$bandwidth,
    0.3
  )

  expect_equal(
    dim(result$score_matrix),
    c(3, length(h_grid))
  )
})


test_that("repeated_cv correctly minimizes an objective", {

  x <- 1:20
  h_grid <- c(0.1, 0.2, 0.3, 0.4)

  mock_score <- function(x, h, fold_id, ...) {
    (h - 0.2)^2
  }

  result <- repeated_cv(
    x = x,
    h_grid = h_grid,
    score_function = mock_score,
    objective = "min",
    k = 5,
    repeats = 4
  )

  expect_equal(
    result$bandwidth,
    0.2
  )

  expect_equal(
    dim(result$score_matrix),
    c(4, length(h_grid))
  )
})


test_that("repeated_cv keeps folds fixed within each repetition", {

  x <- 1:20
  h_grid <- c(0.1, 0.2, 0.3)

  mock_score <- function(x, h, fold_id, ...) {

    # Depends only on the fold assignment, not on h.
    sum(seq_along(fold_id) * fold_id)
  }

  set.seed(1)

  result <- repeated_cv(
    x = x,
    h_grid = h_grid,
    score_function = mock_score,
    objective = "max",
    k = 5,
    repeats = 4
  )

  # Within each repetition every bandwidth should obtain the
  # same value because exactly the same fold_id was used.
  expect_true(
    all(
      apply(
        result$score_matrix,
        1,
        function(z) length(unique(z)) == 1
      )
    )
  )
})