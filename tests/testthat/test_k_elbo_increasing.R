context("elbo")

test_that("elbo is increasing with exponential prior option", {
  set.seed(100)
  sd = 10
  n = 20
  p = 100
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  suppressWarnings(fit.ebmr <- ebmr(X,y,tol = 1e-3,maxiter = 100))

  # Check non-decreasing with tolerance.
  expect_true(all(diff(fit.ebmr$elbo) >= -1e8))
})
