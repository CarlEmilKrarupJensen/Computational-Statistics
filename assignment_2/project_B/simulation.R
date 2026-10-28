# First crossing time for each row of increments; Inf denotes no ruin by horizon.
ruin_times <- function(increments, initial_capital = 30) {
  increments <- as.matrix(increments)
  times <- rep(Inf, nrow(increments))
  for (i in seq_len(nrow(increments))) {
    crossing <- which(initial_capital + cumsum(increments[i, ]) <= 0)
    if (length(crossing) > 0) times[i] <- crossing[1]
  }
  times
}

# Simulate complete paths; return a matrix including time zero for plotting.
simulate_paths <- function(n, theta = 0, initial_capital = 30, horizon = 100,
                           lower = -1.9, upper = 2) {
  validate_model(initial_capital, horizon, lower, upper)
  stopifnot(n >= 1, n == as.integer(n))
  increments <- matrix(r_tilted(n * horizon, theta, lower, upper), nrow = n)
  paths <- matrix(initial_capital, nrow = n, ncol = horizon + 1)
  for (k in seq_len(horizon)) paths[, k + 1] <- paths[, k] + increments[, k]
  list(paths = paths, ruin_times = ruin_times(increments, initial_capital))
}
