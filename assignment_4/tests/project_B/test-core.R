# Run the setup chunk in dose_response_simulation.qmd first, then:
# testthat::test_dir("tests/project_B")
#
# test_that("description", { ... }) groups one check.
# expect_equal(actual, expected, tolerance=...) checks numerical agreement.
# A small tolerance allows for rounding; exact equality is often unsuitable.

library(testthat)

test_that("the gradient and Hessian agree with numerical derivatives", {
  problem <- simulate_doses(100, "continuous", seed = 12)$problem
  theta <- c(0.4, 0.8, 0.1, 1.4)
  h <- 1e-5
  numerical_gradient <- numeric(4)
  numerical_hessian <- matrix(0, 4, 4)

  # Change just one parameter at a time, keeping the other three fixed.
  # The slope between theta-h and theta+h approximates the derivative.
  for (j in 1:4) {
    plus <- theta
    minus <- theta
    plus[j] <- theta[j] + h
    minus[j] <- theta[j] - h

    numerical_gradient[j] <-
      (objective(plus, problem) - objective(minus, problem)) / (2 * h)

    # Differentiating the gradient gives one column of the Hessian.
    numerical_hessian[, j] <-
      (gradient(plus, problem) - gradient(minus, problem)) / (2 * h)
  }

  expect_equal(gradient(theta, problem), numerical_gradient, tolerance = 1e-7)
  # unname removes parameter labels so only the numerical entries are compared.
  expect_equal(unname(hessian(theta, problem)), numerical_hessian, tolerance = 1e-7)
})

test_that("grouping repeated doses leaves the calculations unchanged", {
  # Unequal counts check that grouping uses the correct observation weights.
  x <- exp(c(1, 1, 2, 3, 3, 3))
  y <- c(1.1, 0.9, 0.7, 0.4, 0.3, 0.5)
  problem <- make_problem(x, y)
  grouped <- group_doses(problem)
  theta <- c(2, 1, 0.2, 1.2)

  expect_equal(objective(theta, grouped), objective(theta, problem), tolerance = 1e-10)
  expect_equal(gradient(theta, grouped), gradient(theta, problem), tolerance = 1e-10)
  expect_equal(hessian(theta, grouped), hessian(theta, problem), tolerance = 1e-10)
})

test_that("the full gradient is the mean of the individual gradients", {
  problem <- simulate_doses(20, "continuous", seed = 12)$problem
  theta <- c(0.4, 0.8, 0.1, 1.4)
  total <- numeric(4)

  # Evaluate every individual gradient at the SAME theta, then average.
  # This checks the 1/N scaling. It is not a claim of unbiasedness conditional
  # on earlier updates within a shuffled epoch, where theta changes.
  for (i in 1:problem$N) {
    total <- total + gradient(theta, problem, index = i)
  }
  average <- total / problem$N
  expect_equal(gradient(theta, problem), average, tolerance = 1e-10)
})

test_that("SGD follows the course's shuffled epoch loop", {
  problem <- simulate_doses(20, "continuous", seed = 12)$problem
  start <- initial_parameters(problem)
  theta <- start

  # A short, direct implementation is our reference: shuffle once, then
  # update once per observation. The second epoch uses a smaller rate.
  set.seed(72)
  for (epoch in 1:2) {
    order <- sample.int(problem$N)
    rate <- learning_rate(epoch - 1, gamma0 = 0.01, half_after = 1)
    for (i in order) {
      theta <- theta - rate * gradient(theta, problem, index = i)
    }
  }

  # The same seed gives the same permutations in our solver.
  fit <- stochastic_fit(problem, "sgd", start, batch_size = 1,
                         epochs = 2, gamma0 = 0.01, half_epochs = 1, seed = 72)
  expect_equal(fit$theta, theta, tolerance = 1e-10)
  expect_equal(tail(fit$history$update, 1), 2 * problem$N)
})

test_that("full-gradient methods agree with independent nonlinear least squares", {
  problem <- simulate_doses(500, "continuous", seed = 42)$problem
  reference <- nls(
    y ~ gamma + (delta - gamma) * plogis(alpha - beta * z),
    data = data.frame(y = problem$y, z = problem$z),
    start = as.list(initial_parameters(problem))
  )
  for (method in c("gd", "polyak", "nesterov")) {
    fit <- deterministic_fit(problem, method, maxit = 5000)
    # The safeguards require every accepted full-data step to be non-increasing.
    expect_true(all(diff(fit$history$objective) <= 1e-12))
    expect_equal(objective(fit$theta, problem), mean(residuals(reference)^2),
                 tolerance = 1e-7)
  }
})

test_that("momentum follows the course updates when no safeguard is needed", {
  problem <- simulate_doses(30, "continuous", seed = 12)$problem
  start <- c(0.1, 1, 0.2, 1.2)
  rate <- 0.001
  memory <- 0.9
  for (method in c("polyak", "nesterov")) {
    theta <- start
    previous <- start
    rho <- numeric(4)
    for (n in 1:2) {
      old <- theta
      if (method == "polyak") {
        rho <- memory * rho + (1 - memory) * gradient(theta, problem)
        theta <- theta - rate * rho
      } else {
        point <- theta + memory * (theta - previous)
        theta <- point - rate * gradient(point, problem)
      }
      previous <- old
    }
    fit <- deterministic_fit(problem, method, start, maxit = 2,
                              initial_step = rate, momentum = memory)
    expect_equal(fit$theta, theta, tolerance = 1e-10)
    expect_equal(fit$restarts, 0)
  }
})

test_that("full-batch scheduled GD uses the entire mean gradient", {
  problem <- simulate_doses(20, "continuous", seed = 12)$problem
  start <- initial_parameters(problem)
  fit <- stochastic_fit(problem, "sgd", start, batch_size = problem$N,
                         epochs = 1, gamma0 = 0.01)
  expect_equal(fit$theta, start - 0.01 * gradient(start, problem), tolerance = 1e-10)
})
