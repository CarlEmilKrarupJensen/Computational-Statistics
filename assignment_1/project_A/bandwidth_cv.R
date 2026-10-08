# ============================================================
# Cross-validation utilities
# ============================================================

#' Create approximately balanced random K-fold assignments
make_folds <- function(x, k = 5) {

  n <- length(x)

  if (k < 2 || k > n) {
    stop("k must satisfy 2 <= k <= length(x).")
  }

  sample(
    rep(seq_len(k), length.out = n)
  )

}

#' K-fold likelihood cross-validation score
#'
#' Higher values indicate better predictive log likelihood.
cv_loglik_score <- function(
    x,
    h,
    fold_id,
    kernel = epanechnikov_kernel,
    density_floor = 1e-8
) {

  k <- max(fold_id)

  score <- 0

  for (j in seq_len(k)) {

    test <- x[fold_id == j]
    train <- x[fold_id != j]

    density <- kde_at(
      eval_points = test,
      x = train,
      h = h,
      kernel = kernel
    )

    # Compact-support kernels may assign zero density to held-out
    # observations. The floor prevents log(0), but it introduces
    # a tuning parameter whose sensitivity should be investigated.
    density <- pmax(
      density,
      density_floor
    )

    score <- score + sum(log(density))
  }

  score

}


# ============================================================
# Least-squares cross-validation
# ============================================================

#' Least-squares cross-validation score
#'
#' Computes
#'
#'   LSCV(h) =
#'     integral f_hat_h(x)^2 dx
#'     - (2 / n) * sum_i f_hat_{-i,h}(X_i)
#'
#' where f_hat_{-i,h} is the leave-one-out KDE.
#'
#' Lower values indicate better bandwidths.
#'
#' @param x Numeric vector containing the observed sample.
#' @param h Positive bandwidth.
#' @param kernel Kernel function.
#' @param grid_size Number of points used for numerical integration.
#'
#' @return Numeric LSCV score.
lscv_score <- function(
    x,
    h,
    kernel = epanechnikov_kernel,
    grid_size = 1024
) {

  n <- length(x)

  if (n < 2) {
    stop("LSCV requires at least two observations.")
  }

  if (!is.numeric(h) || length(h) != 1 ||
      !is.finite(h) || h <= 0) {
    stop("h must be a finite positive scalar.")
  }

  # ----------------------------------------------------------
  # First term: integral of the squared KDE
  # ----------------------------------------------------------

  # For the Epanechnikov kernel the estimated density is zero
  # outside [min(x) - h, max(x) + h], so this grid covers the
  # complete support of the KDE.
  grid <- seq(
    min(x) - h,
    max(x) + h,
    length.out = grid_size
  )

  density_grid <- kde_at(
    eval_points = grid,
    x = x,
    h = h,
    kernel = kernel
  )

  # Trapezoidal approximation to integral f_hat_h(x)^2 dx.
  squared_density <- density_grid^2

  integral_term <- sum(
    diff(grid) *
      (head(squared_density, -1) +
         tail(squared_density, -1)) / 2
  )

  # ----------------------------------------------------------
  # Second term: leave-one-out density estimates at X_i
  # ----------------------------------------------------------

  # Evaluate the full KDE at each observed data point.
  density_at_data <- kde_at(
    eval_points = x,
    x = x,
    h = h,
    kernel = kernel
  )

  # The full KDE contains the observation's own contribution
  #
  #     K(0) / (n h).
  #
  # Remove this contribution and renormalize from n to n - 1
  # rather than recomputing n separate KDEs.
  loo_density <- (
    n * density_at_data -
      kernel(0) / h
  ) / (n - 1)

  loo_term <- (2 / n) * sum(loo_density)

  integral_term - loo_term
}

#' K-fold least-squares cross-validation score
#'
#' For each fold, estimates the integrated squared-error criterion
#' using the KDE fitted on the training observations and evaluated
#' on the held-out observations.
#'
#' Lower values indicate better bandwidths.
cv_ls_score <- function(
    x,
    h,
    fold_id,
    kernel = epanechnikov_kernel,
    grid_size = 1024
) {

  k <- max(fold_id)

  fold_scores <- numeric(k)

  for (j in seq_len(k)) {

    test <- x[fold_id == j]
    train <- x[fold_id != j]

    # For a compact-support kernel, the KDE is zero outside this range.
    grid <- seq(
      min(train) - h,
      max(train) + h,
      length.out = grid_size
    )

    density_grid <- kde_at(
      eval_points = grid,
      x = train,
      h = h,
      kernel = kernel
    )

    # First term: integral of the squared training KDE.
    squared_density <- density_grid^2

    integral_term <- sum(
      diff(grid) *
        (
          head(squared_density, -1) +
            tail(squared_density, -1)
        ) / 2
    )

    # Second term: predictive density at held-out observations.
    density_test <- kde_at(
      eval_points = test,
      x = train,
      h = h,
      kernel = kernel
    )

    validation_term <- 2 * mean(density_test)

    fold_scores[j] <- integral_term - validation_term
  }

  mean(fold_scores)
}


#' Repeated K-fold cross-validation
#'
#' Repeats K-fold cross-validation using new random fold assignments
#' for each repetition and averages the resulting objective values.
#'
#' @param x Numeric sample.
#' @param h_grid Candidate bandwidths.
#' @param score_function Cross-validation score function. It must
#'   accept arguments x, h, and fold_id.
#' @param objective Either "max" or "min".
#' @param k Number of folds.
#' @param repeats Number of repeated K-fold partitions.
#' @param ... Additional arguments passed to score_function.
#'
#' @return List containing the selected bandwidth, mean scores,
#'   individual repetition scores, and bandwidth grid.
repeated_cv <- function(
    x,
    h_grid,
    score_function,
    objective = c("max", "min"),
    k = 5,
    repeats = 10,
    ...
) {

  objective <- match.arg(objective)

  score_matrix <- matrix(
    NA_real_,
    nrow = repeats,
    ncol = length(h_grid)
  )

  for (m in seq_len(repeats)) {

    # One fold assignment is generated per repetition and kept
    # fixed across all candidate bandwidths.
    fold_id <- make_folds(x, k)

    for (i in seq_along(h_grid)) {

      score_matrix[m, i] <- score_function(
        x = x,
        h = h_grid[i],
        fold_id = fold_id,
        ...
      )
    }
  }

  mean_score <- colMeans(score_matrix)

  best_index <- switch(
    objective,
    max = which.max(mean_score),
    min = which.min(mean_score)
  )

  list(
    bandwidth = h_grid[best_index],
    mean_score = mean_score,
    score_matrix = score_matrix,
    h_grid = h_grid
  )
}


repeated_cv_loglik <- function(
    x,
    h_grid,
    k = 5,
    repeats = 10,
    kernel = epanechnikov_kernel,
    density_floor = 1e-8
) {

  repeated_cv(
    x = x,
    h_grid = h_grid,
    score_function = cv_loglik_score,
    objective = "max",
    k = k,
    repeats = repeats,
    kernel = kernel,
    density_floor = density_floor
  )
}


repeated_cv_ls <- function(
    x,
    h_grid,
    k = 5,
    repeats = 10,
    kernel = epanechnikov_kernel,
    grid_size = 1024
) {

  repeated_cv(
    x = x,
    h_grid = h_grid,
    score_function = cv_ls_score,
    objective = "min",
    k = k,
    repeats = repeats,
    kernel = kernel,
    grid_size = grid_size
  )
}

cv_loglik_bandwidth <- function(
    x,
    h_grid,
    k = 5,
    kernel = epanechnikov_kernel,
    density_floor = 1e-8
) {

  fold_id <- make_folds(x, k)

  select_bandwidth_grid(
    h_grid = h_grid,
    score_function = function(h) {
      cv_loglik_score(
        x = x,
        h = h,
        fold_id = fold_id,
        kernel = kernel,
        density_floor = density_floor
      )
    },
    objective = "max"
  )
}


lscv_bandwidth <- function(
    x,
    h_grid,
    kernel = epanechnikov_kernel,
    grid_size = 1024
) {

  select_bandwidth_grid(
    h_grid = h_grid,
    score_function = function(h) {
      lscv_score(
        x = x,
        h = h,
        kernel = kernel,
        grid_size = grid_size
      )
    },
    objective = "min"
  )
}