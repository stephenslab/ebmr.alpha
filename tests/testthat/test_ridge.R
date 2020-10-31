context("ridge")

test_that("grr results do not change after simple updates",{
  #grr should not change with further updates as they are optimal for the grr
  set.seed(100)
  sd = 10
  n = 20
  p = 10
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  fit.ebmr = ebmr.init(X,y,1)
  fit.grr = ebmr.update.grr(fit.ebmr, compute_Sigma_full = TRUE, tol=1-15)
  fit.grr = ebmr.update.elbo(fit.grr)

  fit.ebmr = ebmr(X,y,tol=1e-10,maxiter =1000)

  fit.rr = ebmr.update.Sigma.woodbury(fit.grr, compute_Sigma_full = TRUE) # should not move as it should be optimal
  expect_equal(fit.grr,fit.rr,tol=1e-8)

  fit.rr2 = ebmr.update.ebnv.ridge(fit.rr)
  expect_equal(fit.grr,fit.rr2,tol=1e-8)
})


test_that("ridge regression results match simple em",{
  set.seed(100)
  sd = 10
  n = 20
  p = 10
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  fit = ebmr.init(X,y)
  fit.ebmr = ebmr(X,y,tol=1e-10,maxiter =1000)
  fit.em = ridge_em1(y,X,1,1,1000)
  expect_equal(fit.ebmr$residual_variance,fit.em$s2, tol=1e-3)
  expect_equal(fit.ebmr$g$w,fit.em$sb2/fit.em$s2, tol=1e-3)
  expect_true(all(diff(fit.ebmr$elbo)>=-1e8)) #check non-decreasing with tolerance
})

test_that("ridge regression results with admm match simple em",{
  set.seed(100)
  sd = 10
  n = 20
  p = 10
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  suppressWarnings(fit.ebmr <- ebmr(X,y,tol = 1e-12,admm = TRUE))
  fit.em = ridge_em1(y,X,1,1,200)

  expect_equal(fit.ebmr$residual_variance,fit.em$s2,scale = 1,tol = 1e-3)
  expect_equal(fit.ebmr$g$w,fit.em$sb2/fit.em$s2,scale = 1,tol = 1e-3)

  # Check non-decreasing with tolerance.
  expect_true(all(diff(fit.ebmr$elbo) >= -1e8))
})
