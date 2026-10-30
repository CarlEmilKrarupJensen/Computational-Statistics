test_that("expected information matches direct integration of score products", {
 for(nu in c(2,5,30)) {
 v <- 2
 location <- integrate(function(z) ((nu+1)*z/(nu+z^2)/sqrt(v))^2*dt(z,nu),-Inf,Inf)$value
 scale <- integrate(function(z) ((-1+(nu+1)*z^2/(nu+z^2))/(2*v))^2*dt(z,nu),-Inf,Inf)$value
 expect_equal(diag(expected_information(v,nu)),c(location,scale),tolerance=1e-7)
 }
 expect_equal(diag(expected_information(2,1e9)),c(.5,1/8),tolerance=1e-7)
})
test_that("observed information and its transformation agree away from optimum", {
 x <- c(-3,-1,0,1,4,6); mu <- .4; v <- 1.8
 for(p in c("log","original")) expect_equal(observed_information(x,mu,v,5,parameterization=p),
 observed_information(x,mu,v,5,numerical=TRUE,parameterization=p),tolerance=1e-5)
 direct <- -numDeriv::hessian(function(p) t_loglik(x,p[1],p[2],5),c(mu,v))
 expect_equal(observed_information(x,mu,v,5),direct,tolerance=1e-5)
})

test_that("Monte Carlo standard errors use the sigma delta-method transformation", {
  study <- run_monte_carlo(repetitions=2,n=80,nu_values=5,seed=103)
  expect_equal(nrow(study$fits),12L) # two datasets, two starts, three optimizers
  expect_equal(nrow(study$uncertainty),8L) # two datasets, two parameters, two information choices
  expected <- subset(study$uncertainty,information=="Expected")
  location <- subset(expected,parameter=="mu")
  scale <- subset(expected,parameter=="sigma")
  expect_equal(location$se,scale$estimate*sqrt(8/(80*6)),tolerance=1e-10)
  expect_equal(scale$se,scale$estimate*sqrt(8/(2*80*5)),tolerance=1e-10)
  summary <- uncertainty_summary(study)
  expect_equal(summary$empirical_SD[summary$parameter=="mu"],rep(sd(location$estimate),2))
  expect_true(all(summary$valid_SE==2))
  # Zero failures still have a nonzero upper bound for the failure probability.
  timing <- runtime_summary(study)
  expect_true(all(timing$failure_upper>0))
  expect_equal(timing$failure_rate,timing$failures/timing$fits)
  again <- run_monte_carlo(repetitions=2,n=80,nu_values=5,seed=103)
  expect_equal(study$uncertainty,again$uncertainty,tolerance=1e-12)
})
