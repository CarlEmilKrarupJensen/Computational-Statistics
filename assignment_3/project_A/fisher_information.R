# Expected information for N observations, in (mu,v) or (mu,log(v)).
expected_information <- function(v, nu, n=1, parameterization=c("original","log")) {
  parameterization <- match.arg(parameterization); validate_model(c(0,1),0,v,nu)
  if (length(n)!=1 || !is.finite(n) || n<=0) stop("n must be positive")
  result <- diag(c((nu+1)/((nu+3)*v),nu/(2*(nu+3)*v^2)))*n
  if(parameterization=="log") result <- diag(c(1,v)) %*% result %*% diag(c(1,v))
  result
}
# Exact Hessian transformation includes a score term away from a stationary point.
observed_information <- function(x, mu, v, nu, numerical=FALSE, parameterization=c("original","log")) {
  parameterization <- match.arg(parameterization); eta <- c(mu,log(v))
  derivatives <- t_derivatives(x,eta,nu)
  h <- if(numerical) numDeriv::hessian(function(e) t_loglik(x,e[1],exp(e[2]),nu),eta) else derivatives$hessian
  if(parameterization=="original") {
    jacobian <- diag(c(1,1/v)); h <- jacobian %*% h %*% jacobian
    h[2,2] <- h[2,2]-derivatives$gradient[2]/v^2
  }
  -h
}
