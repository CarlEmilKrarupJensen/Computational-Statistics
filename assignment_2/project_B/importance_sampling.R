# Log likelihood f/g for a path with increment sum and specified number of steps.
log_likelihood_ratio <- function(increment_sum, steps, theta,
                                 lower = -1.9, upper = 2) {
  steps * log_mgf(theta, lower, upper) - theta * increment_sum
}

# Summarize independent contributions; distinguish sample and estimator variance.
summarize_contributions <- function(contributions, ruined, steps) {
  n <- length(contributions)
  estimate <- mean(contributions)
  variance <- var(contributions)
  se <- sqrt(variance / n)
  data.frame(n = n, estimate = estimate, standard_error = se,
             ci_lower = estimate - qnorm(0.975) * se,
             ci_upper = estimate + qnorm(0.975) * se,
             contribution_variance = variance, estimator_variance = variance / n,
             relative_standard_error = if (estimate > 0) se / estimate else NA_real_,
             ruin_frequency = mean(ruined), mean_steps = mean(steps))
}

# Estimate ruin with crude MC, full-path IS, or IS stopped at first ruin.
# Optional diagnostics return contributions, actual ruin times, and simulated steps.
# Seeds belong to the caller so independent pilot/final streams remain explicit.
estimate_ruin <- function(n, theta = 0, method = c("stopped", "full", "crude"),
                          initial_capital = 30, horizon = 100,
                          lower = -1.9, upper = 2, diagnostics = FALSE) {
  method <- match.arg(method)
  validate_model(initial_capital, horizon, lower, upper)
  stopifnot(length(n) == 1L, is.finite(n), n >= 2, n == as.integer(n),
            length(theta) == 1L, is.finite(theta))
  if (method == "crude") theta <- 0
  capital <- rep(initial_capital, n)
  times <- rep(Inf, n)
  steps <- integer(n)
  # Full paths continue after ruin; stopped paths draw only for surviving paths.
  # The same time-wise vectorization is used for all methods for comparability.
  for (k in seq_len(horizon)) {
    active <- if (method == "stopped") which(is.infinite(times)) else seq_len(n)
    if (length(active) == 0) break
    capital[active] <- capital[active] + r_tilted(length(active), theta, lower, upper)
    steps[active] <- k
    newly_ruined <- active[is.infinite(times[active]) & capital[active] <= 0]
    times[newly_ruined] <- k
  }
  ruined <- is.finite(times)
  contributions <- numeric(n)
  # Non-ruined paths contribute zero: do not exponentiate unnecessary weights.
  contributions[ruined] <- exp(log_likelihood_ratio(
    capital[ruined] - initial_capital, steps[ruined], theta, lower, upper))
  result <- list(summary = summarize_contributions(contributions, ruined, steps),
                 method = method, theta = theta)
  if (diagnostics) {
    result$contributions <- contributions
    result$ruin_times <- times
    result$steps <- steps
  }
  result
}
