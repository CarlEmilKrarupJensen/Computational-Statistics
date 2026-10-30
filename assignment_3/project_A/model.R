# Shared validation: x is a finite numeric sample; v is scale squared, not variance.
validate_model <- function(x, mu, v, nu) {
  if (!is.numeric(x) || length(x) < 2 || any(!is.finite(x))) stop("x must contain at least two finite numbers")
  if (length(mu) != 1 || !is.finite(mu) || length(v) != 1 || !is.finite(v) || v <= 0 ||
      length(nu) != 1 || !is.finite(nu) || nu <= 0) stop("Require finite mu, v > 0 and nu > 0")
}
# Marginal sample log-likelihood, including all normalizing constants.
t_loglik <- function(x, mu, v, nu) {
  validate_model(x, mu, v, nu)
  sum(dt((x - mu) / sqrt(v), df = nu, log = TRUE)) - length(x) * log(v) / 2
}
# Joint log-density of paired observations under W ~ chi-square(nu).
joint_logdensity <- function(x, w, mu, v, nu) {
  validate_model(x, mu, v, nu)
  if (length(w) != length(x) || any(!is.finite(w) | w <= 0)) stop("Invalid latent observations")
  dchisq(w, nu, log = TRUE) + dnorm(x, mu, sqrt(nu * v / w), log = TRUE)
}
# Gradient and Hessian of log-likelihood in eta = (mu, log(v)).
# q = r^2 / (nu*v + r^2) keeps the derivative expressions compact.
t_derivatives <- function(x, eta, nu) {
  mu <- eta[1]; v <- exp(eta[2]); validate_model(x, mu, v, nu)
  r <- x - mu; a <- nu * v; d <- a + r^2; q <- r^2 / d
  gradient <- c(sum((nu + 1) * r / d), sum(-0.5 + (nu + 1) * q / 2))
  cross <- sum(-(nu + 1) * a * r / d^2)
  hessian <- matrix(c(sum((nu + 1) * (r^2 - a) / d^2), cross,
                      cross, sum(-(nu + 1) * q * (1 - q) / 2)), 2)
  list(gradient = gradient, hessian = hessian)
}
