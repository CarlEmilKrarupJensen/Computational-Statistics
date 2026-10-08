# ============================================================
# Simulation utilities
# ============================================================


#' Integrated squared error
#'
#' Approximates
#'
#'   ISE = integral (f_hat(x) - f(x))^2 dx
#'
#' using trapezoidal integration.
#'
#' @param grid Evaluation grid.
#' @param estimated_density Estimated density values on grid.
#' @param true_density Function evaluating the true density.
#'
#' @return Approximate integrated squared error.
integrated_squared_error <- function(
    grid,
    estimated_density,
    true_density
) {

  true_values <- true_density(grid)

  squared_error <- (
    estimated_density - true_values
  )^2

  sum(
    diff(grid) *
      (
        head(squared_error, -1) +
          tail(squared_error, -1)
      ) / 2
  )
}

#' Evaluate KDE performance for a given bandwidth
#'
#' Fits a KDE using a supplied bandwidth and computes its
#' integrated squared error relative to a known true density.
#'
#' @param x Observed sample.
#' @param h Bandwidth.
#' @param true_density True density function.
#' @param grid Evaluation grid.
#' @param kernel Kernel function.
#'
#' @return Integrated squared error.
evaluate_bandwidth <- function(
    x,
    h,
    true_density,
    grid,
    kernel = epanechnikov_kernel
) {

  estimated_density <- kde_at(
    eval_points = grid,
    x = x,
    h = h,
    kernel = kernel
  )

  integrated_squared_error(
    grid = grid,
    estimated_density = estimated_density,
    true_density = true_density
  )
}


#' Oracle bandwidth on a candidate grid
#'
#' Selects the bandwidth minimizing the true integrated squared
#' error. This is available only in simulation because the true
#' density is known.
#'
#' @param x Observed sample.
#' @param h_grid Candidate bandwidths.
#' @param true_density True density function.
#' @param grid Density evaluation grid.
#' @param kernel Kernel function.
#'
#' @return List containing oracle bandwidth and ISE values.
oracle_bandwidth <- function(
    x,
    h_grid,
    true_density,
    grid,
    kernel = epanechnikov_kernel
) {

  ise <- vapply(
    h_grid,
    function(h) {
      evaluate_bandwidth(
        x = x,
        h = h,
        true_density = true_density,
        grid = grid,
        kernel = kernel
      )
    },
    numeric(1)
  )

  best <- which.min(ise)

  list(
    bandwidth = h_grid[best],
    ise = ise,
    h_grid = h_grid
  )
}


#' Simulation scenarios used in the KDE study
simulation_scenarios <- function() {

  list(

    normal = list(
      label = "Standard normal",

      generate = function(n) {
        rnorm(n)
      },

      density = function(x) {
        dnorm(x)
      },

      range = c(-4, 4)
    ),

    gamma = list(
      label = "Gamma(2, 1)",

      generate = function(n) {
        rgamma(n, shape = 2, rate = 1)
      },

      density = function(x) {
        dgamma(x, shape = 2, rate = 1)
      },

      range = c(0, 10)
    ),

    mixture = list(
      label = "Gaussian mixture",

      generate = function(n) {

        component <- rbinom(
          n,
          size = 1,
          prob = 0.5
        )

        ifelse(
          component == 0,
          rnorm(n, mean = -2, sd = 0.7),
          rnorm(n, mean = 2, sd = 1)
        )
      },

      density = function(x) {

        0.5 * dnorm(
          x,
          mean = -2,
          sd = 0.7
        ) +
          0.5 * dnorm(
            x,
            mean = 2,
            sd = 1
          )
      },

      range = c(-6, 6)
    )

  )
}


#' Apply bandwidth selectors to one dataset
#'
#' @param x Observed sample.
#' @param h_grid Candidate bandwidths for CV methods.
#' @param k Number of folds for likelihood CV.
#' @param repeats Number of repetitions for repeated CV.
#'
#' @return Named vector of selected bandwidths.
select_simulation_bandwidths <- function(
    x,
    h_grid,
    k = 5,
    repeats = 5
) {

  h_amise <- amise_bandwidth(x)

  loglik_result <- repeated_cv_loglik(
    x = x,
    h_grid = h_grid,
    k = k,
    repeats = repeats
  )

  lscv_result <- lscv_bandwidth(
    x = x,
    h_grid = h_grid
  )

  c(
    AMISE = h_amise,
    LogLik_CV = loglik_result$bandwidth,
    LSCV = lscv_result$bandwidth
  )
}


#' Evaluate selected bandwidths against the known density
evaluate_methods <- function(
    x,
    bandwidths,
    true_density,
    grid
) {

  ise <- vapply(
    bandwidths,
    function(h) {

      evaluate_bandwidth(
        x = x,
        h = h,
        true_density = true_density,
        grid = grid
      )

    },
    numeric(1)
  )

  data.frame(
    method = names(bandwidths),
    bandwidth = unname(bandwidths),
    ise = unname(ise)
  )
}


#' Run one simulation replicate
run_simulation_once <- function(
    scenario,
    n,
    h_grid,
    grid_size = 2048,
    k = 5,
    repeats = 5
) {

  x <- scenario$generate(n)

  grid <- seq(
    scenario$range[1],
    scenario$range[2],
    length.out = grid_size
  )

  bandwidths <- select_simulation_bandwidths(
    x = x,
    h_grid = h_grid,
    k = k,
    repeats = repeats
  )

  oracle <- oracle_bandwidth(
    x = x,
    h_grid = h_grid,
    true_density = scenario$density,
    grid = grid
  )

  bandwidths <- c(
    bandwidths,
    Oracle = oracle$bandwidth
  )

  evaluate_methods(
    x = x,
    bandwidths = bandwidths,
    true_density = scenario$density,
    grid = grid
  )
}

#' Repeated KDE simulation study
run_simulation <- function(
    scenario,
    n,
    h_grid,
    B = 30,
    grid_size = 2048,
    k = 5,
    repeats = 5
) {

  results <- vector(
    "list",
    B
  )

  for (b in seq_len(B)) {

    result <- run_simulation_once(
      scenario = scenario,
      n = n,
      h_grid = h_grid,
      grid_size = grid_size,
      k = k,
      repeats = repeats
    )

    result$replicate <- b
    result$n <- n

    results[[b]] <- result
  }

  do.call(
    rbind,
    results
  )
}

#| label: single-example-helper

plot_single_example <- function(
    scenario,
    n = 200,
    h_grid,
    seed = 1
) {

  set.seed(seed)

  x <- scenario$generate(n)

  grid <- seq(
    scenario$range[1],
    scenario$range[2],
    length.out = 1000
  )

  bandwidths <- select_simulation_bandwidths(
    x = x,
    h_grid = h_grid
  )

  oracle <- oracle_bandwidth(
    x = x,
    h_grid = h_grid,
    true_density = scenario$density,
    grid = grid
  )

  bandwidths <- c(
    bandwidths,
    Oracle = oracle$bandwidth
  )

  density_data <- do.call(
    rbind,
    lapply(
      names(bandwidths),
      function(method) {

        data.frame(
          x = grid,
          density = kde_at(
            eval_points = grid,
            x = x,
            h = bandwidths[[method]]
          ),
          method = method
        )
      }
    )
  )

  true_data <- data.frame(
    x = grid,
    density = scenario$density(grid)
  )

  results <- evaluate_methods(
    x = x,
    bandwidths = bandwidths,
    true_density = scenario$density,
    grid = grid
  )

  p <- ggplot() +
    geom_histogram(
      data = data.frame(x = x),
      aes(
        x = x,
        after_stat(density)
      ),
      bins = 30,
      alpha = 0.15
    ) +
    geom_line(
      data = density_data,
      aes(
        x = x,
        y = density,
        color = method
      ),
      linewidth = 0.9
    ) +
    geom_line(
      data = true_data,
      aes(
        x = x,
        y = density
      ),
      linewidth = 0.8,
      linetype = 2,
    ) + 
    labs(
      title = scenario$label,
      subtitle = paste0("n = ", n),
      x = "x",
      y = "Density",
      color = "Method"
    ) + scale_color_brewer(palette = "Dark2") +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom"
    )

  list(
    plot = p,
    sample = x,
    bandwidths = bandwidths,
    results = results
  )
}


compare_bandwidth_selectors_once <- function(
    scenario,
    n = 200,
    h_grid = seq(0.05, 2, length.out = 80),
    grid_size = 2048,
    k = 5,
    repeats = 3
) {

  # Simulate data
  x <- scenario$generate(n)

  # Evaluation grid
  grid <- seq(
    scenario$range[1],
    scenario$range[2],
    length.out = grid_size
  )

  true_density <- scenario$density(grid)

  # -------------------------
  # Your bandwidth selectors
  # -------------------------

  h_amise <- amise_bandwidth(x)

  h_loglik <- repeated_cv_loglik(
    x = x,
    h_grid = h_grid,
    k = k,
    repeats = repeats
  )$bandwidth

  h_lscv <- lscv_bandwidth(
    x = x,
    h_grid = h_grid
  )$bandwidth

  # -------------------------
  # R bandwidth selectors
  # -------------------------

  # Convert Gaussian-SD bandwidth to the support-scale
  # convention of our Epanechnikov kernel.
  h_nrd0 <- sqrt(5) * bw.nrd0(x)

  h_ucv <- sqrt(5) * bw.ucv(x)

  h_sj <- sqrt(5) * bw.SJ(x)

  # -------------------------
  # Oracle bandwidth
  # -------------------------

  oracle <- oracle_bandwidth(
    x = x,
    h_grid = h_grid,
    true_density = scenario$density,
    grid = grid
  )

  h_oracle <- oracle$bandwidth

  bandwidths <- c(
    AMISE = h_amise,
    LogLik_CV = h_loglik,
    LSCV = h_lscv,
    R_nrd0 = h_nrd0,
    R_UCV = h_ucv,
    R_SJ = h_sj,
    Oracle = h_oracle
  )

  # Evaluate every bandwidth using EXACTLY THE SAME KDE
  results <- do.call(
    rbind,
    lapply(
      names(bandwidths),
      function(method) {

        h <- bandwidths[[method]]

        estimate <- kde_at(
          eval_points = grid,
          x = x,
          h = h
        )

        data.frame(
          method = method,
          bandwidth = h,
          ise = integrated_squared_error(
            grid = grid,
            estimated_density = estimate,
            true_density = true_density
          )
        )
      }
    )
  )

  results$oracle_bandwidth <- h_oracle
  results$bandwidth_ratio <-
    results$bandwidth / h_oracle

  oracle_ise <- results$ise[
    results$method == "Oracle"
  ]

  results$ise_ratio <-
    results$ise / oracle_ise

  results
}


run_selector_comparison <- function(
    scenario,
    B = 50,
    n = 200,
    h_grid = seq(0.05, 2, length.out = 80),
    k = 5,
    repeats = 3
) {

  out <- vector(
    "list",
    B
  )

  for (b in seq_len(B)) {

    result <- compare_bandwidth_selectors_once(
      scenario = scenario,
      n = n,
      h_grid = h_grid,
      k = k,
      repeats = repeats
    )

    result$replicate <- b

    out[[b]] <- result
  }

  do.call(
    rbind,
    out
  )
}