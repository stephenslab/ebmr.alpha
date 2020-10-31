test_that("elbo computation matches log likelihood for ridge",{
  set.seed(100)
  sd = 10
  n = 20
  p = 10
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  fit.ebmr = ebmr(X,y,tol=1e-10,maxiter =1000, admm=FALSE)
  expect_equal(elbo(fit.ebmr),compute.ridge.loglik(fit.ebmr))
})
