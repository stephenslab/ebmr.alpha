test_that("elbo computation matches log likelihood for ridge",{
  set.seed(100)
  sd = 10
  n = 20
  p = 30
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  fit.init = ebmr.init(X,y)
  fit.grr = ebmr.update.grr.svd(fit.init,tol=1e-10,maxiter =1000)
  expect_equal(elbo(fit.grr),compute.ridge.loglik(fit.grr))
})
