# Validate the model inputs; return invisibly on success.
validate_model <- function(initial_capital, horizon, lower, upper) {
  stopifnot(length(initial_capital) == 1L, is.finite(initial_capital),
            initial_capital > 0, length(horizon) == 1L,
            is.finite(horizon), horizon >= 1, horizon == as.integer(horizon),
            length(lower) == 1L, length(upper) == 1L,
            is.finite(lower), is.finite(upper), lower < upper)
  invisible(TRUE)
}

# Original uniform density at x, for the supplied support endpoints.
uniform_density <- function(x, lower = -1.9, upper = 2) {
  dunif(x, lower, upper)
}

# Log M(theta), where M is the original increment moment generating function.
# Centering the support gives M(theta) = exp(theta*m) sinh(z)/z.
log_mgf <- function(theta, lower = -1.9, upper = 2) {
  stopifnot(all(is.finite(theta)), lower < upper)
  z <- theta * (upper - lower) / 2
  a <- abs(z)
  correction <- numeric(length(theta))
  small <- a < 1e-3
  correction[small] <- z[small]^2 / 6 - z[small]^4 / 180 + z[small]^6 / 2835
  correction[!small] <- a[!small] + log1p(-exp(-2 * a[!small])) - log(2 * a[!small])
  theta * (lower + upper) / 2 + correction
}

# Log normalizer phi(theta) of the tilted density (not the original MGF).
log_phi <- function(theta, lower = -1.9, upper = 2) {
  log(upper - lower) + log_mgf(theta, lower, upper)
}

# Normalized tilted density at x; theta=0 returns the original uniform density.
tilted_density <- function(x, theta, lower = -1.9, upper = 2) {
  stopifnot(length(theta) == 1L, is.finite(theta))
  result <- exp(theta * x - log_phi(theta, lower, upper))
  result[x < lower | x > upper] <- 0
  result
}

# Draw n tilted increments by inverse transform, with a uniform limit at zero.
r_tilted <- function(n, theta = 0, lower = -1.9, upper = 2) {
  stopifnot(length(theta) == 1L, is.finite(theta), lower < upper)
  u <- runif(n)
  if (theta == 0) return(lower + (upper - lower) * u)
  # Reflect positive tilts to avoid exponentiating a large positive number.
  if (theta > 0) {
    return(upper - log1p(u * expm1(-theta * (upper - lower))) / (-theta))
  }
  lower + log1p(u * expm1(theta * (upper - lower))) / theta
}
