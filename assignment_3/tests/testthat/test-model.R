test_that("joint and marginal densities match their defining distributions", {
 x <- c(-2,.2,3); w <- c(.5,2,7)
 expect_equal(joint_logdensity(x,w,1,2,5),dchisq(w,5,log=TRUE)+dnorm(x,1,sqrt(10/w),log=TRUE))
 expect_equal(t_loglik(x,1,2,5),sum(dt((x-1)/sqrt(2),5,log=TRUE)-log(2)/2))
 marginal <- integrate(function(w) exp(joint_logdensity(rep(.7,length(w)),w,1,2,5)),0,Inf)$value
 expect_equal(marginal,dt((.7-1)/sqrt(2),5)/sqrt(2),tolerance=1e-7)
})
test_that("complete-data analytical solution maximizes the likelihood", {
 x <- c(-2,0,1,4); w <- c(.5,2,3,1)
 analytical <- complete_data_mle(x,w,5)
 numerical <- optim(c(0,0),function(e) -complete_data_loglik(x,w,e[1],exp(e[2]),5),method="BFGS",control=list(reltol=1e-13))
 expect_equal(unname(analytical),c(numerical$par[1],exp(numerical$par[2])),tolerance=1e-5)
})
