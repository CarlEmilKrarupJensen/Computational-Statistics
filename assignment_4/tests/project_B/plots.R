# Plot the actual objective. A dashed reference line avoids assuming the
# best-found numerical solution is a certified lower bound.
plot_convergence <- function(fits, reference) {
  traces <- data.frame()
  for (name in names(fits)) {
    fit <- fits[[name]]
    rows <- data.frame(run = name, method = fit$method,
                       seconds = fit$history$elapsed, mse = fit$history$objective)
    traces <- rbind(traces, rows)
  }
  ggplot2::ggplot(traces, ggplot2::aes(seconds, mse, color = method, group = run)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = reference$value, linetype = "dashed") +
    ggplot2::scale_y_log10() + ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::labs(x = "Elapsed seconds (monitoring included)", y = "Mean squared error",
                 subtitle = "Dashed line: best-found BFGS reference; not a global certificate")
}

# Plot against log-dose so the transition and plateaus are both visible.
plot_curves <- function(dataset, fits) {
  problem <- dataset$problem
  z <- seq(min(problem$z), max(problem$z), length.out = 300)
  truth <- data.frame(z = z, response = dose_response(dataset$truth, z))
  curves <- data.frame()
  for (name in names(fits)) {
    curves <- rbind(curves, data.frame(z = z, response = dose_response(fits[[name]]$theta, z),
                                      method = name))
  }
  # Limit plotted points only; every observation is used for fitting.
  index <- unique(round(seq(1, problem$N, length.out = min(1500, problem$N))))
  points <- data.frame(z = problem$z[index], response = problem$y[index])
  ggplot2::ggplot(points, ggplot2::aes(z, response)) +
    ggplot2::geom_point(alpha = 0.12) +
    ggplot2::geom_line(data = curves, ggplot2::aes(color = method), linewidth = 0.8) +
    ggplot2::geom_line(data = truth, ggplot2::aes(linetype = "Generating curve"),
                        color = "black", linewidth = 0.9) +
    ggplot2::scale_linetype_manual(values = c("Generating curve" = "dashed"), name = NULL) +
    ggplot2::theme_minimal() + ggplot2::labs(x = "log(dose)", y = "Response")
}

# One panel per method avoids hiding overlapping trajectories.
# Color identifies the starting point, consistently across all panels.
plot_slice <- function(problem, fixed, starts, gamma0 = 0.01,
                       batch_size = 1, half_epochs = 5, view = c("contour", "surface")) {
  view <- match.arg(view)
  conditional_fit <- reference_fit(problem, start = fixed, free = 1:2)
  reference <- data.frame(alpha = conditional_fit$theta[1],
                          beta = conditional_fit$theta[2])
  methods <- study_methods()
  labels <- method_labels()
  settings <- study_settings()
  settings$gamma0 <- gamma0
  settings$batch_size <- batch_size
  settings$half_epochs <- half_epochs
  settings$monitor_epochs <- 0.2
  paths <- data.frame()
  endpoints <- data.frame()
  start_labels <- paste0("Start ", seq_len(nrow(starts)), ": (",
                         starts[, 1], ", ", starts[, 2], ")")

  # Keep the same batch size and decay half-time as the main experiment.
  for (j in seq_len(nrow(starts))) {
    start <- fixed
    start[1:2] <- starts[j, ]
    for (method in methods) {
      fit <- study_fit(problem, method, start, settings, seed = 201, free = 1:2)
      h <- fit$history
      trajectory <- data.frame(alpha = h$alpha, beta = h$beta, mse = h$objective,
                                method = method, start = start_labels[j])
      paths <- rbind(paths, trajectory)
      endpoints <- rbind(endpoints, tail(trajectory, 1))
    }
  }
  paths$method <- factor(paths$method, levels = methods)
  endpoints$method <- factor(endpoints$method, levels = methods)
  starting_points <- data.frame(alpha = starts[, 1], beta = starts[, 2],
                                 start = start_labels)

  # All panels use the same axes and exactly the same objective surface.
  ar <- range(c(paths$alpha, reference$alpha), finite = TRUE) + c(-0.5, 0.5)
  br <- range(c(paths$beta, reference$beta), finite = TRUE) + c(-0.1, 0.1)
  surface <- expand.grid(alpha = seq(ar[1], ar[2], length.out = 110),
                         beta = seq(br[1], br[2], length.out = 110))
  surface$mse <- NA_real_
  for (i in seq_len(nrow(surface))) {
    theta <- fixed
    theta[1:2] <- c(surface$alpha[i], surface$beta[i])
    surface$mse[i] <- objective(theta, problem)
  }

  # Optional 3D view of the SAME surface and recorded trajectories.
  # Height is log10(MSE), to make the low-loss valley visible. The solver still
  # minimizes the original MSE; this transform is used only for display.
  if (view == "surface") {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(3, 3), mar = c(1, 1, 3, 1), oma = c(2, 0, 2, 0))
    alpha_grid <- sort(unique(surface$alpha))
    beta_grid <- sort(unique(surface$beta))
    heights <- matrix(log10(pmax(surface$mse, 1e-12)),
                      nrow = length(alpha_grid), ncol = length(beta_grid))
    reference_mse <- objective(conditional_fit$theta, problem)
    # Identical height limits and viewing angle in every panel.
    height_limits <- range(c(heights, log10(pmax(paths$mse, 1e-12)),
                             log10(max(reference_mse, 1e-12))), finite = TRUE)
    start_colors <- rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
                        length.out = length(start_labels))
    for (method in methods) {
      projection <- persp(alpha_grid, beta_grid, heights,
        theta = 40, phi = 30, expand = 0.7, zlim = height_limits,
        col = "#E2EAF0", border = NA, shade = 0.15,
        xlab = "alpha", ylab = "beta", zlab = "log10(MSE)",
        ticktype = "detailed", nticks = 3, cex.axis = 0.65,
        cex.lab = 0.85, main = labels[[method]])
      # trans3d projects each (alpha, beta, log-loss) point onto the figure.
      # Draw trajectories after the surface so they remain visible.
      for (j in seq_along(start_labels)) {
        path <- paths[paths$method == method & paths$start == start_labels[j], ]
        projected <- trans3d(path$alpha, path$beta,
                             log10(pmax(path$mse, 1e-12)), projection)
        lines(projected, col = start_colors[j], lwd = 2)
        points(projected$x[1], projected$y[1], pch = 21,
               bg = "white", col = start_colors[j], cex = 1.2)
        points(tail(projected$x, 1), tail(projected$y, 1),
               pch = 17, col = start_colors[j], cex = 1.1)
      }
      points(trans3d(reference$alpha, reference$beta,
                      log10(max(reference_mse, 1e-12)), projection),
             pch = 8, cex = 1.2)
    }
    mtext(paste("Objective surface; initial stochastic rate =", gamma0),
           outer = TRUE, side = 3, line = 0.3)
    mtext("Circle: start   |   Triangle: last recorded point   |   Star: conditional reference",
           outer = TRUE, side = 1, line = 0.3, cex = 0.8)
    return(invisible(NULL))
  }

  # Log-spaced MSE levels reveal the narrow low-loss valley. This changes only
  # the background color/contour spacing, not the fitted objective or trajectories.
  lower <- max(min(surface$mse), 1e-12)
  upper <- max(max(surface$mse), lower * 1.01)
  levels <- exp(seq(log(lower), log(upper), length.out = 9))
  colors <- rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
                length.out = length(start_labels))
  names(colors) <- start_labels

  ggplot2::ggplot(surface, ggplot2::aes(alpha, beta)) +
    ggplot2::geom_raster(ggplot2::aes(fill = pmax(mse, 1e-12))) +
    ggplot2::scale_fill_gradient(low = "#FAFCFD", high = "#8297A5",
                                 trans = "log10", name = "MSE (log scale)") +
    ggplot2::geom_contour(ggplot2::aes(z = mse), breaks = levels,
                          color = "#667780", linewidth = 0.25, alpha = 0.45) +
    ggplot2::geom_path(data = paths,
      ggplot2::aes(color = start, group = start), linewidth = 0.8,
      arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")) +
    ggplot2::geom_point(data = starting_points, ggplot2::aes(color = start),
                        shape = 21, fill = "white", size = 3, stroke = 1) +
    ggplot2::geom_point(data = endpoints, ggplot2::aes(color = start),
                        shape = 17, size = 2.5) +
    ggplot2::geom_point(data = reference, shape = 8, size = 3, color = "black") +
    ggplot2::scale_color_manual(values = colors, name = NULL) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(
      barwidth = grid::unit(6, "cm"), barheight = grid::unit(0.35, "cm"))) +
    ggplot2::facet_wrap(~method, ncol = 3,
                        labeller = ggplot2::as_labeller(labels)) +
    ggplot2::labs(
      title = paste("2D trajectories: scheduled/AdaGrad initial rate =", gamma0),
      x = expression(alpha), y = expression(beta)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom", legend.box = "vertical")
}

# Compare accuracy against wall time: separate time scales, shared MSE scale.
# Each line is one run; grouped methods are included only when requested.
plot_study <- function(results, include_grouped = FALSE) {
  traces <- data.frame()
  references <- data.frame()
  for (name in names(results)) {
    case <- results[[name]]
    label <- paste(case$dataset$design, "| N =", case$dataset$problem$N)
    references <- rbind(references, data.frame(case = label, mse = case$reference$value))
    for (run in names(case$fits)) {
      if (!include_grouped && startsWith(run, "grouped_")) next
      fit <- case$fits[[run]]
      traces <- rbind(traces, data.frame(case = label, run = run,
        method = fit$method, seconds = fit$history$elapsed,
        mse = fit$history$objective))
    }
  }
  # Preserve the study order: small N on the left, large N on the right.
  panel_order <- unique(references$case)
  traces$case <- factor(traces$case, levels = panel_order)
  references$case <- factor(references$case, levels = panel_order)
  ggplot2::ggplot(traces, ggplot2::aes(seconds, mse, color = method, group = run)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_hline(data = references, ggplot2::aes(yintercept = mse),
                         linetype = "dashed") +
    ggplot2::facet_wrap(~case, ncol = 2, scales = "free_x") + ggplot2::scale_y_log10() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::labs(x = "Elapsed seconds (monitoring included)", y = "Mean squared error",
      subtitle = "Separate time scales; shared MSE scale. Dashed lines: best-found references")
}

# Color gives time, labels give attainment. Grey cells mean nobody reached
# that target within the recorded budget, not a zero-second solution.
plot_accuracy_targets <- function(targets) {
  targets$accuracy <- factor(format(targets$tolerance, scientific = TRUE),
                             levels = format(sort(unique(targets$tolerance), decreasing = TRUE), scientific = TRUE))
  targets$attainment <- paste0(targets$reached, "/", targets$total)
  ggplot2::ggplot(targets, ggplot2::aes(accuracy, method)) +
    ggplot2::geom_tile(ggplot2::aes(fill = pmax(median_seconds, 0.001)), color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = attainment), size = 3) +
    ggplot2::facet_wrap(~case, ncol = 2) +
    ggplot2::scale_fill_gradient(low = "#F0F9E8", high = "#6BAED6",
      trans = "log10", na.value = "grey85", name = "Seconds") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Allowed MSE above the numerical reference (stricter to the right)",
      y = NULL, subtitle = "Labels: runs reaching target / runs attempted; grey: not reached",
      caption = "Color: median first-attainment time among successful runs. Display floor: 1 ms.")
}

# One point per run: lower and further left is better. A common reference is
# used within each dataset, with a display floor for zero/negative differences.
plot_batch_tradeoff <- function(batch_results, accuracy = 1e-4) {
  rows <- batch_results$summary
  rows$batch_label <- factor(rows$batch)
  ggplot2::ggplot(rows, ggplot2::aes(seconds, pmax(difference_from_reference, 1e-8),
                                    color = batch_label, shape = target_reached)) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = accuracy, linetype = "dashed") +
    ggplot2::facet_grid(method ~ case, scales = "free_x") +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16)) +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(x = "Total elapsed seconds (fixed epoch budget)",
      y = "Final MSE minus reference (display floor: 1e-8)", color = "Batch size",
      shape = "Target reached at least once",
      subtitle = "Lower-left is better; dashed line is the accuracy target")
}
