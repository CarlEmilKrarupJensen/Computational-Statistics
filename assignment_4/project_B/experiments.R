# Shared settings: both reports use exactly the same primary experiment.
study_settings <- function() {
  list(sizes = c(300, 15000), epochs = 15, batch_size = 15,
       gamma0 = 0.05, half_epochs = 5, a = 0.75, momentum = 0.9,
       seeds = c(201, 202), monitor_epochs = 0.25,
       maxit = 2000, tolerance = 1e-8, armijo_c = 1e-4, initial_step = 1,
       accuracy = 1e-4)
}

# Prefixes distinguish FULL gradients from SINGLE-observation gradients.
study_methods <- function() {
  c("full_gd", "full_polyak", "full_nesterov", "full_decay",
    "sgd", "polyak", "nesterov", "adagrad")
}
method_labels <- function() {
  c(full_gd = "Full GD: Armijo", full_polyak = "Full Polyak: Armijo",
    full_nesterov = "Full Nesterov: Armijo", full_decay = "Full GD: decay",
    sgd = "SGD", polyak = "Stochastic Polyak",
    nesterov = "Stochastic Nesterov", adagrad = "AdaGrad")
}

study_start <- function(design) {
  if (design == "grid") c(6, 0.8, 0.1, 1.3) else c(-0.5, 0.8, 0.1, 1.3)
}

# Independent numerical reference, never presented as a global certificate.
# free also allows a reference for the conditional two-parameter plots.
reference_fit <- function(problem, start = initial_parameters(problem),
                          free = seq_len(problem$p)) {
  starts <- list(start, start, start, start)
  starts[[2]][1] <- start[1] - 2
  starts[[3]][1] <- start[1] + 2
  starts[[4]][1:2] <- -start[1:2]
  starts[[4]][3:4] <- start[4:3]
  fits <- list()
  rows <- data.frame()
  for (i in seq_along(starts)) {
    # Coefficients outside free remain fixed at the supplied start.
    expand <- function(par) { theta <- start; theta[free] <- par; theta }
    fit <- optim(starts[[i]][free],
      function(par) objective(expand(par), problem),
      function(par) gradient(expand(par), problem)[free],
      method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
    fits[[i]] <- expand(fit$par)
    rows <- rbind(rows, data.frame(start = i, objective = fit$value,
      convergence_code = fit$convergence,
      gradient_norm = sqrt(sum(gradient(fits[[i]], problem)[free]^2))))
  }
  best <- which.min(rows$objective)
  list(theta = fits[[best]], value = rows$objective[best], runs = rows)
}

# The only place where study controls are passed to the solvers.
study_fit <- function(problem, method, start, settings, seed = 201,
                       free = seq_len(problem$p)) {
  if (method %in% c("full_gd", "full_polyak", "full_nesterov")) {
    fit <- deterministic_fit(problem, sub("full_", "", method), start,
      maxit = settings$maxit, tolerance = settings$tolerance,
      armijo_c = settings$armijo_c, initial_step = settings$initial_step,
      momentum = settings$momentum, free = free)
  } else {
    # One batch of N observations is one full-gradient update per epoch.
    # It shares the decay schedule with SGD but is NOT an Armijo method.
    batch <- if (method == "full_decay") problem$N else settings$batch_size
    algorithm <- if (method == "full_decay") "sgd" else method
    fit <- stochastic_fit(problem, algorithm, start, batch_size = batch,
      epochs = settings$epochs, gamma0 = settings$gamma0,
      half_epochs = settings$half_epochs, a = settings$a,
      momentum = settings$momentum, seed = seed,
      monitor_epochs = settings$monitor_epochs, free = free)
  }
  fit$method <- unname(method_labels()[method])
  fit$seed <- if (startsWith(method, "full_")) NA_integer_ else seed
  fit
}

# One large dataset per design; smaller samples are prefixes of it.
study_data <- function(settings) {
  cases <- list()
  for (design in c("grid", "continuous")) {
    full <- simulate_doses(max(settings$sizes), design, seed = 2026)
    for (N in settings$sizes) {
      dataset <- full
      index <- seq_len(N)
      dataset$problem <- make_problem(exp(full$problem$z[index]), full$problem$y[index])
      cases[[paste(design, N, sep = "_")]] <- dataset
    }
  }
  cases
}

# Reuse this runner for size, start, rate and geometry experiments.
# No rate is chosen after inspecting a method's result.
compare_methods <- function(problem, start, settings, methods = study_methods()) {
  fits <- list()
  for (method in methods) {
    seeds <- if (startsWith(method, "full_")) settings$seeds[1] else settings$seeds
    for (seed in seeds) {
      name <- if (startsWith(method, "full_")) method else paste(method, seed, sep = "_")
      fits[[name]] <- study_fit(problem, method, start, settings, seed)
    }
  }
  fits
}

run_study <- function(settings = study_settings()) {
  cases <- study_data(settings)
  results <- list()
  for (name in names(cases)) {
    dataset <- cases[[name]]
    problem <- dataset$problem
    start <- study_start(dataset$design)
    fits <- compare_methods(problem, start, settings)
    reference_problem <- problem
    preprocessing <- 0
    if (dataset$design == "grid") {
      preprocessing <- system.time(grouped <- group_doses(problem))[["elapsed"]]
      reference_problem <- grouped
      for (method in c("full_gd", "full_polyak", "full_nesterov")) {
        fit <- study_fit(grouped, method, start, settings)
        fit$history$elapsed <- fit$history$elapsed + preprocessing
        fit$method <- paste("Grouped", fit$method)
        fits[[paste0("grouped_", method)]] <- fit
      }
    }
    results[[name]] <- list(dataset = dataset, start = start, fits = fits,
      reference = reference_fit(reference_problem), preprocessing = preprocessing)
  }
  results
}

# Safe summaries also retain failed runs rather than hiding them.
fit_summary <- function(fits, problem, reference, accuracy = 1e-4) {
  out <- data.frame()
  for (name in names(fits)) {
    fit <- fits[[name]]
    h <- fit$history
    value <- if (all(is.finite(fit$theta))) objective(fit$theta, problem) else NA_real_
    norm <- if (is.finite(value)) sqrt(sum(gradient(fit$theta, problem)^2)) else NA_real_
    reached <- which(is.finite(h$objective) & h$objective <= reference$value + accuracy)
    out <- rbind(out, data.frame(run = name, method = fit$method,
      seed = if (is.null(fit$seed)) NA_integer_ else fit$seed, status = fit$status,
      seconds = tail(h$elapsed, 1), updates = tail(h$update, 1), objective = value,
      difference_from_reference = value - reference$value, gradient_norm = norm,
      seconds_to_target = if (length(reached)) h$elapsed[reached[1]] else NA_real_,
      restarts = if (is.null(fit$restarts)) NA_integer_ else fit$restarts))
  }
  out
}

study_summary <- function(results, settings) {
  out <- data.frame()
  for (name in names(results)) {
    case <- results[[name]]
    rows <- fit_summary(case$fits, case$dataset$problem, case$reference, settings$accuracy)
    rows$design <- case$dataset$design
    rows$N <- case$dataset$problem$N
    rows$reference_mse <- case$reference$value
    out <- rbind(out, rows)
  }
  rownames(out) <- NULL
  out
}

# Conditioning is LOCAL and depends on parameter units. The exact Hessian can
# be indefinite away from a minimum; do not report a positive-definite condition
# number in that case. 2 J' W J is a positive-semidefinite curvature diagnostic.
curvature_summary <- function(theta, problem) {
  J <- response_gradient(theta, problem$z)
  GN <- 2 * crossprod(J, J * problem$weights)
  e <- eigen(hessian(theta, problem), symmetric = TRUE, only.values = TRUE)$values
  q <- eigen(GN, symmetric = TRUE, only.values = TRUE)$values
  ratio <- function(values) {
    if (min(values) <= max(values) * 1e-10) Inf else max(values) / min(values)
  }
  data.frame(hessian_min = min(e), hessian_max = max(e),
    hessian_condition = if (min(e) > 0) ratio(e) else NA_real_,
    jacobian_curvature_condition = ratio(q))
}

# Change only the true slope: retain doses, noise realization, center and plateaus.
# A shallow sigmoid can look almost linear over the observed range, weakening
# identification. A steep sigmoid concentrates information near its transition.
geometry_data <- function(N = 1500) {
  set.seed(2026)
  z <- rnorm(N, sd = 2)
  noise <- rnorm(N, sd = 0.1)
  cases <- list()
  for (slope in c(0.35, 1.2, 4)) {
    truth <- c(alpha = 0, beta = slope, gamma = 0.2, delta = 1.2)
    cases[[paste0("slope_", slope)]] <- list(
      problem = make_problem(exp(z), dose_response(truth, z) + noise),
      truth = truth, design = "continuous")
  }
  cases
}

# Small focused experiments, rather than a large factorial tuning grid.
# Every returned case has the same structure as the main sample-size study.
run_sensitivity <- function(settings, experiment = "starts") {
  experiment <- match.arg(experiment, c("starts", "stochastic_rate", "armijo_step", "geometry"))
  base <- study_data(settings)[["grid_300"]]
  if (experiment == "geometry") {
    datasets <- geometry_data()
    values <- seq_along(datasets)
  } else if (experiment == "starts") {
    starts <- list(near = c(7.5, 0.9, 0.15, 1.25),
                   baseline = study_start("grid"), saturated = c(20, 2, 0.5, 2.3))
    values <- names(starts)
  } else if (experiment == "stochastic_rate") {
    values <- c(0.01, 0.05, 0.1)
  } else {
    values <- c(0.1, 1, 5)
  }
  results <- list()
  for (value in values) {
    controls <- settings
    dataset <- base
    start <- study_start("grid")
    methods <- study_methods()
    label <- paste(experiment, value, sep = " = ")
    if (experiment == "starts") start <- starts[[value]]
    if (experiment == "stochastic_rate") {
      controls$gamma0 <- value
      methods <- c("full_decay", "sgd", "polyak", "nesterov", "adagrad")
    }
    if (experiment == "armijo_step") {
      controls$initial_step <- value
      methods <- c("full_gd", "full_polyak", "full_nesterov")
    }
    if (experiment == "geometry") {
      dataset <- datasets[[value]]
      start <- study_start("continuous")
      label <- names(datasets)[value]
    }
    reference_problem <- if (dataset$design == "grid") group_doses(dataset$problem) else dataset$problem
    results[[label]] <- list(dataset = dataset, start = start,
      fits = compare_methods(dataset$problem, start, controls, methods),
      reference = reference_fit(reference_problem), preprocessing = 0)
  }
  results
}

# Report actual local geometry, not a presumed ordering of "easy" and "hard".
geometry_summary <- function(results) {
  rows <- data.frame()
  for (name in names(results)) {
    case <- results[[name]]
    for (location in c("truth", "start", "reference")) {
      theta <- switch(location, truth = case$dataset$truth,
                       start = case$start, reference = case$reference$theta)
      rows <- rbind(rows, cbind(data.frame(case = name, location = location,
        gradient_norm = sqrt(sum(gradient(theta, case$dataset$problem)^2))),
        curvature_summary(theta, case$dataset$problem)))
    }
  }
  rows
}

# Batch-size experiment: same data, start, shuffle seeds and epoch schedule.
# Return summaries and histories so the report can examine both final accuracy
# and progress over time. This does not tune a separate rate for each batch size.
run_batch_study <- function(results, settings, batches = c(1, 16, 64, 256)) {
  rows <- data.frame()
  traces <- data.frame()
  for (case_name in names(results)) {
    case <- results[[case_name]]
    problem <- case$dataset$problem
    for (batch in batches) {
      controls <- settings
      controls$batch_size <- batch
      for (method in c("sgd", "polyak", "nesterov", "adagrad")) {
        for (seed in settings$seeds) {
          fit <- study_fit(problem, method, case$start, controls, seed)
          summary <- fit_summary(list(run = fit), problem, case$reference, settings$accuracy)
          summary$case <- case_name
          summary$batch <- batch
          summary$target_reached <- is.finite(summary$seconds_to_target)
          rows <- rbind(rows, summary)
          h <- fit$history
          traces <- rbind(traces, data.frame(case = case_name, method = fit$method,
            seed = seed, batch = batch, seconds = h$elapsed,
            passes = h$passes, mse = h$objective,
            reference_difference = h$objective - case$reference$value))
        }
      }
    }
  }
  list(summary = rows, traces = traces)
}

# Time to several requested accuracies. Include non-attainment explicitly.
# The median is among successful runs only; the reached/total label must be
# shown alongside it so a method is not rewarded for hiding unsuccessful runs.
accuracy_targets <- function(results, tolerances = c(1e-3, 1e-4, 1e-5)) {
  runs <- data.frame()
  for (case_name in names(results)) {
    case <- results[[case_name]]
    for (name in names(case$fits)) {
      if (startsWith(name, "grouped_")) next
      fit <- case$fits[[name]]
      for (tolerance in tolerances) {
        h <- fit$history
        reached <- which(is.finite(h$objective) &
                           h$objective <= case$reference$value + tolerance)
        runs <- rbind(runs, data.frame(case = case_name, method = fit$method,
          tolerance = tolerance,
          seconds = if (length(reached)) h$elapsed[reached[1]] else NA_real_))
      }
    }
  }
  groups <- split(runs, interaction(runs$case, runs$method, runs$tolerance, drop = TRUE))
  rows <- lapply(groups, function(g) {
    good <- is.finite(g$seconds)
    data.frame(case = g$case[1], method = g$method[1], tolerance = g$tolerance[1],
      reached = sum(good), total = nrow(g),
      median_seconds = if (any(good)) median(g$seconds[good]) else NA_real_)
  })
  do.call(rbind, rows)
}

# bench stores all repeated durations in a list-column. Expand it for a
# beeswarm plot. Keep GC iterations: they are part of the measured workload.
benchmark_samples <- function(timings) {
  rows <- data.frame()
  for (i in seq_len(nrow(timings))) {
    rows <- rbind(rows, data.frame(implementation = as.character(timings$expression[i]),
      seconds = as.numeric(timings$time[[i]])))
  }
  rows
}
