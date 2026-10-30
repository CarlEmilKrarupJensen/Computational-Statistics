# BFGS uses optim's own stopping rule, then the shared score check certifies accuracy.
# History records objective evaluations (including line-search trials), not iterations.
fit_bfgs <- function(x, nu, start=c(median(x), log(IQR(x)^2/4)), tol=1e-7, maxit=10000) {
  validate_model(x, start[1], exp(start[2]), nu)
  if (tol <= 0 || maxit < 1) stop("Invalid convergence controls")
  tic <- proc.time()[3]; history <- data.frame(); nf <- 0; ng <- 0
  fn <- function(eta) {
    nf <<- nf+1
    value <- tryCatch(-t_loglik(x, eta[1], exp(eta[2]), nu), error=function(e) Inf)
    history <<- rbind(history, data.frame(iteration=nf, time=unname(proc.time()[3]-tic), loglik=-value, mu=eta[1], v=exp(eta[2])))
    value
  }
  gr <- function(eta) { ng <<- ng+1; -t_derivatives(x, eta, nu)$gradient }
  result <- optim(start, fn, gr, method="BFGS", control=list(reltol=min(1e-14,tol^2), maxit=maxit))
  eta <- result$par
  list(mu=unname(eta[1]), v=unname(exp(eta[2])), loglik=-result$value, iterations=NA_integer_,
       converged=result$convergence==0 && score_error(x,eta,nu)<=tol,
       evaluations=c(objective=nf,gradient=ng,hessian=0), elapsed=unname(proc.time()[3]-tic), history=history)
}
# Damped Newton ascent. Indefinite information uses a gradient ascent direction.
fit_newton <- function(x, nu, start=c(median(x), log(IQR(x)^2/4)), tol=1e-7, maxit=10000) {
  validate_model(x, start[1], exp(start[2]), nu)
  if (tol <= 0 || maxit < 1) stop("Invalid convergence controls")
  tic <- proc.time()[3]; eta <- start; ll <- t_loglik(x,eta[1],exp(eta[2]),nu); nf <- 1; nd <- 0
  history <- data.frame(iteration=0,time=0,loglik=ll,mu=eta[1],v=exp(eta[2]))
  converged <- FALSE; reason <- "maximum iterations"; iteration <- 0
  for (k in seq_len(maxit)) {
    derivatives <- t_derivatives(x,eta,nu); nd <- nd+1; g <- derivatives$gradient
    if (score_error(x,eta,nu)<=tol) { converged <- TRUE; reason <- "score tolerance"; break }
    information <- -derivatives$hessian
    direction <- tryCatch({chol(information); solve(information,g)}, error=function(e) g/length(x))
    if (any(!is.finite(direction)) || sum(direction*g)<=0) { reason <- "unsuitable direction"; break }
    alpha <- 1; accepted <- FALSE
    for (j in 1:50) {
      candidate <- eta+alpha*direction
      candidate_ll <- tryCatch(t_loglik(x,candidate[1],exp(candidate[2]),nu),error=function(e) -Inf)
      nf <- nf+1
      if (is.finite(candidate_ll) && candidate_ll >= ll+1e-4*alpha*sum(g*direction)) { accepted <- TRUE; break }
      alpha <- alpha/2
    }
    if (!accepted) { reason <- "line search failed"; break }
    eta <- candidate; ll <- candidate_ll; iteration <- k
    history <- rbind(history,data.frame(iteration=k,time=unname(proc.time()[3]-tic),loglik=ll,mu=eta[1],v=exp(eta[2])))
  }
  converged <- score_error(x,eta,nu)<=tol
  list(mu=unname(eta[1]),v=unname(exp(eta[2])),loglik=ll,iterations=iteration,converged=converged,
       reason=if(converged) "score tolerance" else reason,evaluations=c(objective=nf,gradient=nd,hessian=nd),
       elapsed=unname(proc.time()[3]-tic),history=history)
}
