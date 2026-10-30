# Complete-data log-likelihood for observed pairs (x,w).
complete_data_loglik <- function(x, w, mu, v, nu) sum(joint_logdensity(x, w, mu, v, nu))
# Analytical maximizer. Returns location and scale squared.
complete_data_mle <- function(x, w, nu) {
  joint_logdensity(x, w, 0, 1, nu) # validate the paired sample
  mu <- sum(w * x) / sum(w)
  v <- sum(w * (x - mu)^2) / (length(x) * nu)
  if (v <= 0) stop("No positive-scale maximizer for a constant sample")
  c(mu = mu, v = v)
}
