# ============================================================
# Kernel density estimation
# ============================================================

#' Evaluate a kernel density estimator
#'
#' Evaluates
#'
#'   f_hat(z) = 1/(n h) sum_i K((z - X_i)/h)
#'
#' at arbitrary evaluation points.
#'
#' @param eval_points Numeric vector of points where the KDE is evaluated.
#' @param x Numeric vector containing the observed sample.
#' @param h Positive bandwidth.
#' @param kernel Kernel function. Defaults to the Epanechnikov kernel.
#'
#' @return Numeric vector of estimated density values.
kde_at <- function(
    eval_points,
    x,
    h,
    kernel = epanechnikov_kernel
) {

  if (!is.numeric(x) || length(x) < 1) {
    stop("x must be a non-empty numeric vector.")
  }

  if (any(!is.finite(x))) {
    stop("x must contain only finite values.")
  }

  if (!is.numeric(h) || length(h) != 1 || !is.finite(h) || h <= 0) {
    stop("h must be a finite positive scalar.")
  }

  vapply(
    eval_points,
    function(z) {

      # Standardized distances from the evaluation point
      # to every observation in the sample.
      u <- (z - x) / h

      sum(kernel(u)) / (length(x) * h)

    },
    numeric(1)
  )

}


#' Kernel density estimate on an automatically generated grid
#'
#' @param x Numeric sample.
#' @param h Positive bandwidth.
#' @param n Number of grid points.
#' @param kernel Kernel function.
#' @param extend Number of bandwidths by which the grid is extended
#'   beyond the observed data range.
#'
#' @return List containing evaluation grid, density values,
#'   and bandwidth.
kernel_density_estimator <- function(
    x,
    h,
    n = 512,
    kernel = epanechnikov_kernel,
    extend = 1
) {

  if (n < 2) {
    stop("n must be at least 2.")
  }

  data_range <- range(x)

  grid <- seq(
    data_range[1] - extend * h,
    data_range[2] + extend * h,
    length.out = n
  )

  density <- kde_at(
    eval_points = grid,
    x = x,
    h = h,
    kernel = kernel
  )

  list(
    x = grid,
    y = density,
    bandwidth = h
  )

}