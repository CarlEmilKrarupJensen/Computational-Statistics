select_bandwidth_grid <- function(
    h_grid,
    score_function,
    objective = c("max", "min")
) {

  objective <- match.arg(objective)

  scores <- vapply(
    h_grid,
    score_function,
    numeric(1)
  )

  best_index <- switch(
    objective,
    max = which.max(scores),
    min = which.min(scores)
  )

  list(
    bandwidth = h_grid[best_index],
    scores = scores
  )

}


