# E[W | x]: W|x is Gamma(shape=(nu+1)/2, rate=(1+r^2/(nu*v))/2).
conditional_weights <- function(x, mu, v, nu) {
  validate_model(x, mu, v, nu)
  (nu + 1) / (1 + ((x - mu) / sqrt(nu * v))^2)
}
# One E-step followed by the analytical M-step.
em_step <- function(x, mu, v, nu) complete_data_mle(x, conditional_weights(x, mu, v, nu), nu)
# Common accuracy measure: mean score in dimensionless location/log-scale coordinates.
score_error <- function(x, eta, nu) {
  g <- t_derivatives(x, eta, nu)$gradient / length(x)
  max(abs(g * c(exp(eta[2] / 2), 1)))
}
# Fit with a common score tolerance and a relative likelihood-change check.
# history includes the starting point; time is elapsed seconds since entry.
fit_em <- function(x, nu, start = c(median(x), log(IQR(x)^2 / 4)), tol = 1e-7,
                   maxit = 10000, keep_history = TRUE, optimized = FALSE) {
  validate_model(x, start[1], exp(start[2]), nu)
  if (tol <= 0 || maxit < 1) stop("Invalid convergence controls")
  tic <- proc.time()[3]; eta <- start; ll <- t_loglik(x, eta[1], exp(eta[2]), nu)
  history <- data.frame(iteration=0, time=0, loglik=ll, mu=eta[1], v=exp(eta[2]))
  converged <- FALSE
  for (iteration in seq_len(maxit)) {
    if (optimized) {
      w <- conditional_weights(x, eta[1], exp(eta[2]), nu)
      mu <- sum(w*x)/sum(w)
      estimate <- c(mu, sum(w*(x-mu)^2)/(length(x)*nu))
    } else estimate <- em_step(x, eta[1], exp(eta[2]), nu)
    candidate <- c(estimate[1], log(estimate[2]))
    new_ll <- t_loglik(x, candidate[1], exp(candidate[2]), nu)
    if (new_ll < ll - 1e-10 * (1 + abs(ll))) stop("EM likelihood decreased")
    change <- abs(new_ll-ll)/(1+abs(ll)); eta <- candidate; ll <- new_ll
    history <- rbind(history, data.frame(iteration, time=unname(proc.time()[3]-tic), loglik=ll, mu=eta[1], v=exp(eta[2])))
    if (score_error(x, eta, nu) <= tol && change <= tol) { converged <- TRUE; break }
  }
  list(mu=unname(eta[1]), v=unname(exp(eta[2])), loglik=unname(ll), iterations=iteration,
       converged=converged, evaluations=c(objective=iteration+1, gradient=iteration, hessian=0),
       elapsed=unname(proc.time()[3]-tic), history=if(keep_history) history else NULL)
}
