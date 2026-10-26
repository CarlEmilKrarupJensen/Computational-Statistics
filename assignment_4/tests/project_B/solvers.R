# Assignment B: the objective is nonconvex. Small gradients mean stationarity,
# not a guarantee of a local or global minimum.
# Lecture schedule: gamma_n = gamma0 * K / (K + n^a).
# n counts completed epochs; half_after is the number before the rate has halved.
learning_rate <- function(n, gamma0, half_after, a = 0.75) {
  K <- half_after^a
  gamma0 * K / (K + n^a)
}

# Random reshuffling, as in the course code:
#   outer loop: shuffle once per epoch;
#   inner loop: visit consecutive mini-batches of that permutation.
# Every observation is used exactly once in each completed epoch.
# Only the update rule differs between the four methods.
stochastic_fit <- function(problem, method = "sgd",
                           start = initial_parameters(problem), batch_size = 32,
                           epochs = 30, gamma0 = 0.1, half_epochs = 10,
                           a = 0.75, momentum = 0.9, epsilon = 1e-8,
                           seed = 1, monitor_epochs = 1,
                           free = seq_len(problem$p)) {
  method <- match.arg(method, c("sgd", "polyak", "nesterov", "adagrad"))
  stopifnot(length(start) == problem$p, all(is.finite(start)),
            batch_size >= 1, batch_size <= problem$N, batch_size %% 1 == 0,
            is.finite(epochs), epochs >= 1, epochs %% 1 == 0,
            is.finite(monitor_epochs), monitor_epochs > 0,
            gamma0 > 0, half_epochs > 0, a > 0.5, a <= 1,
            momentum >= 0, momentum < 1, epsilon > 0,
            length(free) > 0, !anyDuplicated(free), all(free %in% seq_len(problem$p)))
  # Shuffle individual observations, not compressed dose groups.
  stopifnot(!problem$grouped)
  set.seed(seed)

  N <- problem$N
  batch_starts <- seq(1, N, by = batch_size)
  batches_per_epoch <- length(batch_starts)

  # monitor_epochs affects recording only, not sampling or the learning rate.
  # Values below 1 allow intermediate points in trajectory plots.
  # Values >= 1 record every ceiling(monitor_epochs) epochs, plus the final epoch.
  record_every <- max(1, ceiling(monitor_epochs * batches_per_epoch))
  records_per_epoch <- floor((batches_per_epoch - 1) / record_every) + 1
  history <- matrix(NA_real_, 1 + epochs * records_per_epoch, 4 + problem$p)
  colnames(history) <- c("update", "passes", "elapsed", "objective",
                         "alpha", "beta", "gamma", "delta")

  theta <- start
  # Initialize ONCE, outside the epoch loop. All memory persists across epochs.
  # Polyak: rho is the exponentially weighted average of gradients.
  # Nesterov: theta_previous stores the previous parameter vector.
  # In the course slides momentum is called beta; here that name would conflict
  # with the dose-response slope, so the argument remains "momentum".
  rho <- numeric(problem$p)
  theta_previous <- start
  squared_gradients <- numeric(problem$p)
  status <- "budget_exhausted"
  update <- 0
  processed <- 0
  row <- 1
  started <- proc.time()[["elapsed"]]
  value <- objective(theta, problem)
  history[row, ] <- c(0, 0, proc.time()[["elapsed"]] - started, value, theta)

  for (epoch in seq_len(epochs)) {
    # A fresh random permutation: every index occurs exactly once.
    permutation <- sample.int(N)

    # Same base rate for every batch in this epoch.
    # Use elapsed epochs (epoch - 1), so epoch 1 starts at gamma0.
    # With half_epochs=10, epoch 11 begins at gamma0/2.
    rate <- learning_rate(epoch - 1, gamma0, half_epochs, a)

    for (batch in seq_along(batch_starts)) {
      first <- batch_starts[batch]
      last <- min(first + batch_size - 1, N)
      index <- permutation[first:last]

      # The last batch may be smaller. gradient() divides by length(index),
      # so it averages over the ACTUAL batch size rather than the nominal one.
      point <- theta
      if (method == "nesterov") {
        # Course notation: y_n = theta_{n-1} + beta_n*(theta_{n-1}-theta_{n-2}).
        # Here point is y_n (not an observed response); beta_n = momentum = 0.9.
        point <- theta + momentum * (theta - theta_previous)
      }
      g <- gradient(point, problem, index)

      # "free" selects the parameters allowed to change in the 2D experiments.
      if (method == "sgd") {
        theta[free] <- theta[free] - rate * g[free]
      }
      if (method == "polyak") {
        # Course convention: rho_n = beta*rho_{n-1} + (1-beta)*g_n.
        # The gradient is evaluated at the current theta, not a look-ahead point.
        rho[free] <- momentum * rho[free] + (1 - momentum) * g[free]
        theta[free] <- theta[free] - rate * rho[free]
      }
      if (method == "nesterov") {
        # Course convention: theta_n = y_n - gamma*g(y_n).
        # Save the old theta BEFORE replacing it, for the next look-ahead.
        theta_previous <- theta
        theta[free] <- point[free] - rate * g[free]
      }
      if (method == "adagrad") {
        # AdaGrad keeps gamma0 fixed but adapts each coordinate after EVERY batch.
        # It does not use the extra epoch-based decay schedule.
        squared_gradients[free] <- squared_gradients[free] + g[free]^2
        theta[free] <- theta[free] - gamma0 * g[free] /
          (sqrt(squared_gradients[free] + epsilon))
      }

      update <- update + 1
      processed <- processed + length(index)
      if (any(!is.finite(theta))) {
        status <- "nonfinite_iterate"
        break
      }

      # Count actual observations processed. At an epoch boundary, passes is
      # exactly the epoch number, even when N is not divisible by batch_size.
      intermediate_record <- monitor_epochs < 1 && batch %% record_every == 0
      epoch_record <- batch == batches_per_epoch &&
        (epoch %% ceiling(monitor_epochs) == 0 || epoch == epochs)
      if (intermediate_record || epoch_record) {
        value <- objective(theta, problem)
        row <- row + 1
        history[row, ] <- c(update, processed / N,
                            proc.time()[["elapsed"]] - started, value, theta)
        if (!is.finite(value)) {
          status <- "nonfinite_objective"
          break
        }
      }
    }
    # A failure must leave BOTH loops, not start another epoch.
    if (status != "budget_exhausted") break
  }

  history <- as.data.frame(history[seq_len(row), , drop = FALSE])
  list(theta = theta, history = history, method = method, status = status)
}

# Armijo search from a supplied base point. Try initial_step, then halve it.
# The direction must point downhill at THIS base point.
armijo_search <- function(theta, value, g, direction, problem, free,
                          initial_step, armijo_c) {
  step <- initial_step
  for (attempt in 1:60) {
    proposal <- theta
    proposal[free] <- theta[free] + step * direction
    proposed_value <- objective(proposal, problem)
    if (is.finite(proposed_value) &&
        proposed_value <= value + armijo_c * step * sum(g * direction)) {
      return(list(accepted = TRUE, theta = proposal, value = proposed_value,
                  step = step, evaluations = attempt))
    }
    step <- step / 2
  }
  list(accepted = FALSE, evaluations = 60)
}

# Full-batch GD, Polyak gradient averaging, or Nesterov look-ahead.
# Momentum plus Armijo is a safeguarded variant, not the fixed-step theorem
# from the convex course slides. Restarting removes unhelpful momentum.
deterministic_fit <- function(problem, method = "gd",
                              start = initial_parameters(problem), maxit = 2000,
                              tolerance = 1e-8, free = seq_len(problem$p),
                              armijo_c = 1e-4, initial_step = 1, momentum = 0.9) {
  method <- match.arg(method, c("gd", "polyak", "nesterov"))
  stopifnot(length(start) == problem$p, all(is.finite(start)), maxit >= 1,
            tolerance > 0, armijo_c > 0, armijo_c < 1,
            is.finite(initial_step), initial_step > 0,
            momentum >= 0, momentum < 1,
            length(free) > 0, !anyDuplicated(free),
            all(free %in% seq_len(problem$p)))
  theta <- start
  previous <- start
  rho <- numeric(length(free))
  status <- "maxit"
  restarts <- 0
  evaluations <- 1
  gradient_evaluations <- 0
  history <- matrix(NA_real_, maxit + 1, 6 + problem$p)
  colnames(history) <- c("update", "elapsed", "objective", "step",
                         "restarts", "objective_evaluations",
                         "alpha", "beta", "gamma", "delta")
  started <- proc.time()[["elapsed"]]
  value <- objective(theta, problem)
  row <- 1
  history[row, ] <- c(0, proc.time()[["elapsed"]] - started, value,
                      0, restarts, evaluations, theta)

  for (n in seq_len(maxit)) {
    g <- gradient(theta, problem)[free]
    gradient_evaluations <- gradient_evaluations + 1
    if (max(abs(g)) < tolerance) {
      status <- "stationary"
      break
    }
    base <- theta
    base_value <- value
    base_gradient <- g
    direction <- -g
    restart <- FALSE

    if (method == "polyak") {
      # Course convention: an exponentially weighted AVERAGE of gradients.
      rho <- momentum * rho + (1 - momentum) * g
      direction <- -rho
      # Past gradients can point uphill at the current parameter vector.
      # Restart the average at zero, then incorporate the current gradient.
      if (sum(g * direction) >= 0) {
        rho <- (1 - momentum) * g
        direction <- -rho
        restarts <- restarts + 1
      }
    }
    if (method == "nesterov") {
      # Extrapolate first, then evaluate the gradient at the look-ahead point.
      base[free] <- theta[free] + momentum * (theta[free] - previous[free])
      base_value <- objective(base, problem)
      base_gradient <- gradient(base, problem)[free]
      evaluations <- evaluations + 1
      gradient_evaluations <- gradient_evaluations + 1
      direction <- -base_gradient
    }

    # Nesterov's search is based at the look-ahead point. GD/Polyak search
    # from the current iterate. The trial step resets on every outer iteration.
    if (is.finite(base_value) && all(is.finite(direction))) {
      trial <- armijo_search(base, base_value, base_gradient, direction,
                              problem, free, initial_step, armijo_c)
      evaluations <- evaluations + trial$evaluations
    } else {
      trial <- list(accepted = FALSE)
    }

    # Armijo at the look-ahead point alone need not decrease H(theta).
    # If the extrapolated step increases it, discard momentum and retry as GD.
    # Also retry a failed momentum line search as GD.
    restart <- method != "gd" &&
      (!trial$accepted || trial$value > value)
    if (restart) {
      restarts <- restarts + 1
      rho[] <- 0
      trial <- armijo_search(theta, value, g, -g, problem, free,
                              initial_step, armijo_c)
      evaluations <- evaluations + trial$evaluations
    }
    if (!trial$accepted) {
      status <- "line_search_failed"
      break
    }
    old_theta <- theta
    theta <- trial$theta
    value <- trial$value
    # After a restart the next Nesterov look-ahead has zero displacement.
    previous <- if (restart) theta else old_theta
    row <- row + 1
    history[row, ] <- c(n, proc.time()[["elapsed"]] - started, value,
                        trial$step, restarts, evaluations, theta)
  }
  if (max(abs(gradient(theta, problem)[free])) < tolerance) status <- "stationary"
  history <- as.data.frame(history[seq_len(row), , drop = FALSE])
  list(theta = theta, history = history, method = paste("full", method),
       status = status, restarts = restarts,
       objective_evaluations = evaluations,
       gradient_evaluations = gradient_evaluations)
}
