test_that("analytical derivatives agree with independent numerical derivatives", {
 x <- c(-3,-1,.1,.4,1,2,7); eta <- c(.3,log(1.7))
 fun <- function(e) t_loglik(x,e[1],exp(e[2]),2)
 d <- t_derivatives(x,eta,2)
 expect_equal(d$gradient,numDeriv::grad(fun,eta),tolerance=1e-6)
 expect_equal(d$hessian,numDeriv::hessian(fun,eta),tolerance=1e-5)
})
test_that("three optimizers agree across fixed degrees of freedom", {
 for(nu in c(2,5,30)) {
 x <- simulate_t_data(500,1,2,nu,seed=42+nu); f <- fit_all(x,nu)
 expect_true(all(sapply(f,`[[`,"converged")))
 expect_lt(diff(range(sapply(f,`[[`,"loglik"))),1e-7)
 expect_lt(diff(range(sapply(f,`[[`,"mu"))),1e-5)
 expect_lt(diff(range(sapply(f,`[[`,"v"))),1e-5)
 }
})
test_that("numerical optimizers report invalid inputs and incomplete runs", {
 x <- c(-2,-.5,0,1,4)
 for(fun in list(fit_newton,fit_bfgs)) {
 expect_error(fun(x,5,start=c(0,Inf)))
 expect_false(fun(x,5,start=c(10,2),maxit=1)$converged)
 }
})

test_that("contamination experiment uses paired data and equivalent starts", {
  study <- run_contamination_study(repetitions=2,n=100,nu_values=c(2,5),fractions=c(0,.05),seed=901)
  expect_equal(nrow(study$fits),48L)
  groups <- split(study$fits,interaction(study$fits[c("replication","fraction","scenario")],drop=TRUE))
  expect_true(all(vapply(groups,function(g) length(unique(g$start_mu))==1,logical(1))))
  expect_true(all(vapply(groups,function(g) length(unique(g$start_logv))==1,logical(1))))
  expect_true(all(study$fits$replacements==100*study$fits$fraction))
  # Reconstruct the exact contaminated sample and verify a reported likelihood.
  x <- simulate_t_data(100,0,1,5,seed=902)
  positions <- sample.int(100); outliers <- rnorm(100,15,1)
  x[positions[1:5]] <- outliers[1:5]
  fit <- subset(study$fits,replication==1 & fraction==.05 & nu==5 & scenario=="Data-based" & method=="EM")
  expect_equal(fit$loglik,t_loglik(x,fit$mu,fit$v,5),tolerance=1e-10)
  summary <- contamination_summary(study)
  expect_equal(sum(summary$fits),48L)
  expect_equal(summary$failure_rate,summary$failures/summary$fits)
  agreement <- contamination_agreement(study)
  expect_equal(nrow(agreement),8L)
  expect_true(all(agreement$certified>=1))
  expect_lt(max(agreement$loglik_spread),1e-7)
})
