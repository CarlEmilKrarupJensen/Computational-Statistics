# Draw from the marginal model; optional seed makes an experiment reproducible.
simulate_t_data <- function(n, mu=0, v=1, nu=5, seed=NULL) {
  validate_model(c(0,1),mu,v,nu)
  if(n<2) stop("n must be at least two")
  if(!is.null(seed)) set.seed(seed)
  mu+sqrt(v)*rt(n,nu)
}
