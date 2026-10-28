# Run independent pilot experiments over theta, methods, and sample sizes.
variance_experiment <- function(theta_grid, sample_sizes = c(5000, 20000), seed = 2101) {
  set.seed(seed)
  rows <- list()
  index <- 0L
  for (n in sample_sizes) {
    for (method in c("full", "stopped")) {
      for (theta in theta_grid) {
        index <- index + 1L
        row <- estimate_ruin(n, theta, method)$summary
        row$theta <- theta
        row$method <- method
        rows[[index]] <- row
      }
    }
  }
  do.call(rbind, rows)
}

# Select the best grid point for each method from a pilot table and objective.
select_theta <- function(results, objective = "contribution_variance") {
  groups <- split(results, results$method)
  do.call(rbind, lapply(groups, function(group) {
    group[which.min(group[[objective]]), , drop = FALSE]
  }))
}

# Benchmark fixed-size batches after one warm-up; keep every repetition.
# RNG initialization and garbage collection occur outside the timed expression.
benchmark_methods <- function(theta_grid, n = 20000, repetitions = 5, seed = 3101,
                              methods = c("full", "stopped")) {
  stopifnot(repetitions >= 2)
  jobs <- expand.grid(theta = theta_grid, method = methods,
                      repetition = seq_len(repetitions), stringsAsFactors = FALSE)
  set.seed(seed)
  jobs <- jobs[sample.int(nrow(jobs)), ] # Random order reduces systematic timing drift.
  estimate_ruin(n, theta_grid[1], methods[1])
  rows <- vector("list", nrow(jobs))
  for (i in seq_len(nrow(jobs))) {
    set.seed(seed + i)
    gc()
    elapsed <- system.time(result <- estimate_ruin(n, jobs$theta[i], jobs$method[i]))[["elapsed"]]
    row <- result$summary
    row$theta <- jobs$theta[i]
    row$method <- jobs$method[i]
    row$repetition <- jobs$repetition[i]
    row$elapsed <- elapsed
    row$cost_per_sample <- elapsed / n
    row$variance_cost <- row$contribution_variance * row$cost_per_sample
    rows[[i]] <- row
  }
  do.call(rbind, rows)
}

# Aggregate repetitions: mean pilot variance and median time, with timing spread.
summarize_benchmarks <- function(results) {
  groups <- split(results, interaction(results$method, results$theta, drop = TRUE))
  do.call(rbind, lapply(groups, function(group) {
    data.frame(method = group$method[1], theta = group$theta[1],
               contribution_variance = mean(group$contribution_variance),
               estimate = mean(group$estimate), mean_steps = mean(group$mean_steps),
               elapsed = median(group$elapsed),
               elapsed_min = min(group$elapsed), elapsed_max = max(group$elapsed),
               cost_per_sample = median(group$cost_per_sample),
               variance_cost = mean(group$contribution_variance) * median(group$cost_per_sample))
  }))
}

# Extrapolate sample size and time for a target RSE from independent pilot estimates.
precision_budget <- function(results, probability, target_rse = 0.05) {
  stopifnot(probability > 0, target_rse > 0)
  results$required_n <- pmax(2, ceiling(results$contribution_variance /
                                        (target_rse * probability)^2))
  results$required_seconds <- results$required_n * results$cost_per_sample
  results
}

# Profile repeated batches for a fixed duration; return Rprof's self-time table.
profile_method <- function(method, theta, n = 20000, seconds = 2,
                            file = tempfile(fileext = ".out"), seed = 4101) {
  set.seed(seed)
  start <- proc.time()[["elapsed"]]
  Rprof(file, interval = 0.01)
  on.exit(Rprof(NULL), add = TRUE)
  repeat {
    estimate_ruin(n, theta, method)
    if (proc.time()[["elapsed"]] - start >= seconds) break
  }
  Rprof(NULL)
  summaryRprof(file)$by.self
}

# Profile fixed batches by source line; return times joined to the actual R code.
# Source modules with keep.source=TRUE before calling, so R retains line locations.
profile_source_lines <- function(method, theta, n = 50000, batches = 50,
                                 seed = 4201, source_dir = "project_B",
                                 file = tempfile(fileext = ".out")) {
  stopifnot(batches >= 1, batches == as.integer(batches))
  if (is.null(attr(estimate_ruin, "srcref"))) {
    stop("Source the R modules with keep.source = TRUE before line profiling.")
  }
  # Warm-up and garbage collection are outside the sampled workload.
  estimate_ruin(n, theta, method)
  gc()
  set.seed(seed)
  Rprof(file, interval = 0.01, line.profiling = TRUE)
  on.exit(Rprof(NULL), add = TRUE)
  for (batch in seq_len(batches)) estimate_ruin(n, theta, method)
  Rprof(NULL)
  profile <- summaryRprof(file, lines = "show")
  rows <- profile$by.line
  locations <- rownames(rows)
  # Report simulation and sampler statements; retain the whole-run denominator.
  keep <- grepl("(importance_sampling|model)[.]R#[0-9]+$", locations)
  rows <- rows[keep, , drop = FALSE]
  locations <- locations[keep]
  if (nrow(rows) == 0) stop("No simulation source lines were sampled.")
  rows$source_file <- basename(sub("#[0-9]+$", "", locations))
  rows$line <- as.integer(sub(".*#", "", locations))
  rows$code <- character(nrow(rows))
  for (i in seq_len(nrow(rows))) {
    source <- readLines(file.path(source_dir, rows$source_file[i]))
    rows$code[i] <- trimws(source[rows$line[i]])
  }
  rows$location <- paste0(rows$source_file, ":", rows$line)
  rows <- rows[order(rows$self.time, decreasing = TRUE), , drop = FALSE]
  rownames(rows) <- NULL
  list(lines = rows, sampling_time = profile$sampling.time,
       sample_interval = profile$sample.interval, batches = batches, n = n)
}
