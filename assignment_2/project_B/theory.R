# Derivative K'(theta): the exact mean increment under the tilted distribution.
tilted_mean <- function(theta, lower = -1.9, upper = 2) {
  half_width <- (upper - lower) / 2
  z <- theta * half_width
  result <- numeric(length(theta))
  small <- abs(z) < 1e-3
  result[small] <- z[small] / 3 - z[small]^3 / 45 + 2 * z[small]^5 / 945
  result[!small] <- 1 / tanh(z[!small]) - 1 / z[!small]
  (lower + upper) / 2 + half_width * result
}

