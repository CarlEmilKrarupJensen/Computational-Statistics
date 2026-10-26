# Parameter order throughout: alpha, beta, gamma, delta.
# Work with z = log(x) to avoid repeatedly evaluating logarithms.
dose_response <- function(theta, z) {
  q <- plogis(theta[1] - theta[2] * z)
  theta[3] + (theta[4] - theta[3]) * q
}

# First derivatives of f with respect to the four parameters (one row per dose).
response_gradient <- function(theta, z) {
  q <- plogis(theta[1] - theta[2] * z)
  w <- q * plogis(-theta[1] + theta[2] * z)
  difference <- theta[4] - theta[3]
  cbind(alpha = difference * w, beta = -difference * z * w,
        gamma = 1 - q, delta = q)
}

# Mean squared error, including the within-dose variation for grouped data.
objective <- function(theta, problem) {
  residual <- dose_response(theta, problem$z) - problem$y
  sum(problem$weights * residual^2) + problem$constant
}

# Full gradient, or the average over the current observation mini-batch.
gradient <- function(theta, problem, index = NULL) {
  z <- problem$z
  y <- problem$y
  weights <- problem$weights
  if (!is.null(index)) {
    # Sampling is from original observations, not uniformly from dose groups.
    stopifnot(!problem$grouped)
    z <- z[index]
    y <- y[index]
    weights <- rep(1 / length(index), length(index))
  }
  residual <- dose_response(theta, z) - y
  J <- response_gradient(theta, z)
  as.vector(2 * crossprod(J, weights * residual))
}

# Exact Hessian of MSE: 2 sum_i w_i [J_i J_i' + (f_i-y_i) Hessian(f_i)].
# The residual term can make this matrix indefinite away from a minimum.
hessian <- function(theta, problem) {
  z <- problem$z
  q <- plogis(theta[1] - theta[2] * z)
  w <- q * plogis(-theta[1] + theta[2] * z)
  t <- w * (1 - 2 * q)
  difference <- theta[4] - theta[3]
  residual <- dose_response(theta, z) - problem$y
  J <- response_gradient(theta, z)
  H <- 2 * crossprod(J, J * problem$weights)

  # Nonzero entries of the response Hessian; fill the symmetric entries too.
  R <- matrix(0, 4, 4)
  a <- 2 * problem$weights * residual
  R[1, 1] <- sum(a * difference * t)
  R[1, 2] <- R[2, 1] <- sum(-a * difference * z * t)
  R[2, 2] <- sum(a * difference * z^2 * t)
  R[1, 3] <- R[3, 1] <- sum(-a * w)
  R[1, 4] <- R[4, 1] <- sum(a * w)
  R[2, 3] <- R[3, 2] <- sum(a * z * w)
  R[2, 4] <- R[4, 2] <- sum(-a * z * w)
  H + R
}
