test_that("E and M steps match a hand-calculated example", {
 # nu=2, mu=0, v=1: weights are 3,2,1 for residuals 0,1,2.
 x <- c(0,1,2); w <- c(3,2,1)
 expect_equal(conditional_weights(x,0,1,2),w)
 expect_equal(unname(em_step(x,0,1,2)),c(2/3,5/9))
})
test_that("EM increases likelihood and streamlined implementation agrees", {
 x <- simulate_t_data(500,1,2,5,seed=42)
 a <- fit_em(x,5); b <- fit_em(x,5,optimized=TRUE)
 expect_true(a$converged); expect_true(all(diff(a$history$loglik)>=-1e-9))
 expect_equal(c(a$mu,a$v,a$loglik),c(b$mu,b$v,b$loglik),tolerance=1e-12)
 expect_false(fit_em(x,5,start=c(8,0),maxit=1)$converged)
 expect_error(fit_em(x,5,start=c(0,-Inf)))
 expect_error(fit_em(rep(1,10),5,start=c(1,0)))
})
