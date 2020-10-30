context("elbo")

test_that("elbo is increasing with k option", {
  set.seed(100)
  sd = 10
  n = 200
  p = 1000
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  fit.ebmr = ebmr(X,y,tol = 1e-10,maxiter = 100,admm = TRUE,k = 10)

  # Check non-decreasing with tolerance.
  expect_true(all(diff(fit.ebmr$elbo) >= -1e8)) 
})
