# ============================================================
# AMISE plug-in bandwidth
# ============================================================

#' Robust scale estimate used in Silverman's rule
#'
#' Computes
#'
#'   min(sd(x), IQR(x) / 1.34).
silverman_scale <- function(x) {

  min(
    sd(x),
    IQR(x) / 1.34
  )

}

#' Silverman rule-of-thumb bandwidth
#'
#' Used here as a pilot bandwidth.
silverman_bandwidth <- function(x) {

  n <- length(x)

  sigma_hat <- silverman_scale(x)

  0.9 * sigma_hat * n^(-1 / 5)

}


#' Second derivative of the standard Gaussian density
gaussian_second_derivative <- function(u) {

  (u^2 - 1) * dnorm(u)

}


#' Numerical estimate of R(f'')
#'
#' Computes the integrated squared second derivative of a
#' Gaussian pilot density using numerical integration.
#'
#' This implementation is intentionally computationally expensive
#' and is retained primarily as a reference implementation.
Rf2_gaussian_numerical <- function(x, r) {

  n <- length(x)

  pair_integral <- function(i, j) {

    integrand <- function(z) {

      gaussian_second_derivative(
        (z - x[i]) / r
      ) *
        gaussian_second_derivative(
          (z - x[j]) / r
        )

    }

    integrate(
      integrand,
      lower = -Inf,
      upper = Inf
    )$value

  }

  indices <- expand.grid(
    i = seq_len(n),
    j = seq_len(n)
  )

  values <- mapply(
    pair_integral,
    indices$i,
    indices$j
  )

  # Each second derivative contributes a factor r^{-3};
  # integration over z contributes one factor r.
  sum(values) / (n^2 * r^6)

}

#' Closed-form Gaussian convolution term
#'
#' Appears when integrating the product of two second derivatives
#' of Gaussian kernels.
gaussian_second_derivative_convolution <- function(d) {

  (d^4 - 12 * d^2 + 12) *
    exp(-d^2 / 4) /
    (32 * sqrt(pi))

}


#' Closed-form estimate of R(f'')
#'
#' Uses pairwise standardized distances to avoid performing
#' numerical integration for every observation pair.
Rf2_gaussian_closed <- function(x, r) {

  n <- length(x)

  # Matrix of pairwise standardized distances.
  # This vectorizes the O(n^2) pairwise calculation, but also
  # requires O(n^2) memory.
  D <- outer(x, x, "-") / r

  C <- gaussian_second_derivative_convolution(D)

  sum(C) / (n^2 * r^5)

}


#' AMISE plug-in bandwidth
#'
#' Estimates the asymptotically optimal KDE bandwidth
#'
#'   h = [R(K) / {n mu_2(K)^2 R(f'')}]^(1/5).
#'
#' A Gaussian pilot density is used to estimate R(f'').
#'
#' @param x Numeric sample.
#' @param pilot_method Either "closed" or "numerical".
#'
#' @return Estimated AMISE-optimal bandwidth.
amise_bandwidth <- function(
    x,
    pilot_method = c("closed", "numerical")
) {

  pilot_method <- match.arg(pilot_method)

  n <- length(x)

  # Pilot bandwidth used to estimate the unknown roughness R(f'').
  r <- silverman_bandwidth(x)

  Rf2 <- switch(
    pilot_method,
    closed = Rf2_gaussian_closed(x, r),
    numerical = Rf2_gaussian_numerical(x, r)
  )

  RK <- RK_epanechnikov()
  mu2 <- mu2_epanechnikov()

  (
    RK /
      (n * mu2^2 * Rf2)
  )^(1 / 5)

}