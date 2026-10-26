# Store original observations. The weights produce a mean, not a sum of squares.
make_problem <- function(x, y) {
  stopifnot(length(x) == length(y), length(x) > 0, all(is.finite(x)),
            all(x > 0), all(is.finite(y)))
  list(z = log(x), y = y, weights = rep(1 / length(y), length(y)),
       constant = 0, N = length(y), p = 4, grouped = FALSE)
}

# Simulate the two designs described in the assignment.
# Grid doses: exp(1), ..., exp(15); the transition is centered at log(x)=8.
# Continuous doses: log(x) ~ N(0, omega^2); transition centered at log(x)=0.
simulate_doses <- function(N = 3000, design = "grid", seed = 2026,
                           sigma = 0.1, omega = 2) {
  stopifnot(N >= 15, sigma > 0, omega > 0)
  design <- match.arg(design, c("grid", "continuous"))
  set.seed(seed)
  if (design == "grid") {
    z <- rep(1:15, length.out = N)
    theta <- c(alpha = 8, beta = 1, gamma = 0.2, delta = 1.2)
  } else {
    z <- rnorm(N, sd = omega)
    theta <- c(alpha = 0, beta = 1.2, gamma = 0.2, delta = 1.2)
  }
  y <- dose_response(theta, z) + rnorm(N, sd = sigma)
  list(problem = make_problem(exp(z), y), truth = theta, design = design,
       sigma = sigma, omega = omega)
}

# Exact compression of repeated doses. For each dose retain its count, mean,
# and sum of squared deviations from that mean. No rounding or binning is used.
group_doses <- function(problem) {
  stopifnot(!problem$grouped)
  doses <- sort(unique(problem$z))
  group <- match(problem$z, doses)
  counts <- tabulate(group, nbins = length(doses))
  sums <- as.vector(rowsum(problem$y, group, reorder = TRUE))
  means <- sums / counts
  within <- sum((problem$y - means[group])^2)
  list(z = doses, y = means, weights = counts / problem$N,
       constant = within / problem$N, N = problem$N, p = 4, grouped = TRUE)
}

# A data-based decreasing-curve start, independent of the true parameters.
# Quantiles provide distinct initial plateaus; alpha/beta locates the transition.
initial_parameters <- function(problem) {
  lower <- as.numeric(quantile(problem$y, 0.1))
  upper <- as.numeric(quantile(problem$y, 0.9))
  if (upper == lower) upper <- lower + 0.1
  c(alpha = median(problem$z), beta = 1, gamma = lower, delta = upper)
}
