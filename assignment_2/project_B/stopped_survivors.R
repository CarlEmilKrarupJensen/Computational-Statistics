# Estimate ruin by stopping paths and maintaining their surviving indices.
# Same inputs and output as estimate_ruin(method="stopped"), without method selection.
# Only active-index bookkeeping differs; weighting and summaries are unchanged.
estimate_ruin_survivors <- function(n, theta = 0, initial_capital = 30,
                                    horizon = 100, lower = -1.9, upper = 2,
                                    diagnostics = FALSE) {
  validate_model(initial_capital, horizon, lower, upper)
  stopifnot(length(n) == 1L, is.finite(n), n >= 2, n == as.integer(n),
            length(theta) == 1L, is.finite(theta))
  capital <- rep(initial_capital, n)
  times <- rep(Inf, n)
  steps <- integer(n)
  active <- seq_len(n)

  for (k in seq_len(horizon)) {
    if (length(active) == 0) break
    capital[active] <- capital[active] + r_tilted(length(active), theta, lower, upper)
    steps[active] <- k
    # Every active path is a survivor, so no is.infinite() check is needed.
    hit <- capital[active] <= 0
    times[active[hit]] <- k
    # Keep the survivors in their original order. This also preserves the
    # assignment of future random draws, enabling exact paired comparisons.
    active <- active[!hit]
  }

  ruined <- is.finite(times)
  contributions <- numeric(n)
  contributions[ruined] <- exp(log_likelihood_ratio(
    capital[ruined] - initial_capital, steps[ruined], theta, lower, upper))
  result <- list(summary = summarize_contributions(contributions, ruined, steps),
                 method = "stopped", theta = theta)
  if (diagnostics) {
    result$contributions <- contributions
    result$ruin_times <- times
    result$steps <- steps
  }
  result
}

# Benchmark original and survivor-index stopped estimators on identical random draws.
# Return one row per timed batch; warm-up, GC, seeds, and checks are not timed.
benchmark_stopped_versions <- function(theta_grid = c(-0.55, -0.30, -0.10, 0),
                                       n = 150000, repetitions = 15, seed = 5201) {
  stopifnot(repetitions >= 2)
  implementations <- c("Original", "Survivor indices")
  for (theta in theta_grid) {
    estimate_ruin(n, theta, "stopped")
    estimate_ruin_survivors(n, theta)
  }
  set.seed(seed)
  jobs <- expand.grid(theta = theta_grid, repetition = seq_len(repetitions))
  jobs <- jobs[sample.int(nrow(jobs)), ]
  # Randomize both pair order and within-pair order before resetting RNG seeds.
  original_first <- sample(c(TRUE, FALSE), nrow(jobs), replace = TRUE)
  rows <- list()
  index <- 0L
  for (i in seq_len(nrow(jobs))) {
    order <- if (original_first[i]) implementations else rev(implementations)
    pair_results <- list()
    for (implementation in order) {
      gc()
      set.seed(seed + i)
      if (implementation == "Original") {
        elapsed <- system.time(result <- estimate_ruin(n, jobs$theta[i], "stopped"))[["elapsed"]]
      } else {
        elapsed <- system.time(result <- estimate_ruin_survivors(n, jobs$theta[i]))[["elapsed"]]
      }
      pair_results[[implementation]] <- result$summary
      index <- index + 1L
      row <- result$summary
      row$theta <- jobs$theta[i]
      row$repetition <- jobs$repetition[i]
      row$pair <- i
      row$implementation <- implementation
      row$elapsed <- elapsed
      row$cost_per_sample <- elapsed / n
      row$variance_cost <- row$contribution_variance * row$cost_per_sample
      rows[[index]] <- row
    }
    # Exact equality is stronger than agreement within Monte Carlo error.
    stopifnot(identical(pair_results[[1]], pair_results[[2]]))
  }
  do.call(rbind, rows)
}

# Summarize paired runtimes; speedup > 1 means survivor indices were faster.
summarize_stopped_benchmark <- function(results) {
  original <- subset(results, implementation == "Original")
  improved <- subset(results, implementation == "Survivor indices")
  pairs <- merge(original, improved, by = c("theta", "pair", "repetition"),
                 suffixes = c("_original", "_improved"))
  groups <- split(pairs, pairs$theta)
  do.call(rbind, lapply(groups, function(group) {
    data.frame(theta = group$theta[1], n = group$n_original[1],
               repetitions = nrow(group),
               original_seconds = median(group$elapsed_original),
               improved_seconds = median(group$elapsed_improved),
               median_paired_speedup = median(group$elapsed_original / group$elapsed_improved),
               median_paired_time_saved_pct = 100 * median(1 - group$elapsed_improved / group$elapsed_original),
               mean_steps = mean(group$mean_steps_original),
               summaries_identical = all(group$estimate_original == group$estimate_improved &
                 group$contribution_variance_original == group$contribution_variance_improved))
  }))
}

# Benchmark candidate tilts using only the survivor-index implementation.
# Return per-batch variances and costs, compatible with summarize_benchmarks().
benchmark_survivor_tilts <- function(theta_grid, n = 150000, repetitions = 9,
                                     seed = 7201) {
  stopifnot(length(theta_grid) > 0, all(is.finite(theta_grid)), repetitions >= 2)
  theta_grid <- sort(unique(round(theta_grid, 8)))
  estimate_ruin_survivors(n, theta_grid[1])
  jobs <- expand.grid(theta = theta_grid, repetition = seq_len(repetitions))
  set.seed(seed)
  jobs <- jobs[sample.int(nrow(jobs)), ]
  rows <- vector("list", nrow(jobs))
  for (i in seq_len(nrow(jobs))) {
    gc()
    set.seed(seed + i)
    elapsed <- system.time(result <- estimate_ruin_survivors(n, jobs$theta[i]))[["elapsed"]]
    row <- result$summary
    row$method <- "stopped"
    row$theta <- jobs$theta[i]
    row$repetition <- jobs$repetition[i]
    row$elapsed <- elapsed
    row$cost_per_sample <- elapsed / n
    row$variance_cost <- row$contribution_variance * row$cost_per_sample
    rows[[i]] <- row
  }
  do.call(rbind, rows)
}

# Normalize variance and time to a preselected reference tilt; return relative cost.
# R = (v*c)/(v_ref*c_ref), so R<1 predicts less time for equal target precision.
relative_survivor_efficiency <- function(results, reference_theta) {
  reference <- results[abs(results$theta - reference_theta) < 1e-7, , drop = FALSE]
  if (nrow(reference) != 1) stop("Exactly one reference row is required.")
  stopifnot(reference$contribution_variance > 0, reference$cost_per_sample > 0)
  results$relative_variance <- results$contribution_variance / reference$contribution_variance
  results$relative_runtime <- results$cost_per_sample / reference$cost_per_sample
  results$relative_cost <- results$relative_variance * results$relative_runtime
  results$predicted_time_saved_pct <- 100 * (1 - results$relative_cost)
  results[order(results$theta), , drop = FALSE]
}
