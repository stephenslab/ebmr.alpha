test_that("ridge regression results match simple em", {
  set.seed(100)
  sd = 10
  n = 20
  p = 10
  X = matrix(rnorm(n*p),ncol=p)
  btrue = rnorm(p)
  y = X %*% btrue + sd*rnorm(n)
  fit.ebmr = ebmr(X,y,tol=1e-10,maxiter =1000, admm=FALSE)
  fit.em = ridge_em1(y,X,1,1,1000)
  expect_equal(fit.ebmr$residual_variance,fit.em$s2, tol=1e-3)
  expect_equal(fit.ebmr$g,fit.em$sb2/fit.em$s2, tol=1e-3)
  expect_true(all(diff(fit.ebmr$elbo)>=-1e8)) #check non-decreasing with tolerance
})

test_that("ridge regression results with admm match simple em", {
  set.seed(100)
  sd = 10
  n = 20
  p = 10
  X = matrix(rnorm(n*p),ncol=p)
  btrue = rnorm(p)
  y = X %*% btrue + sd*rnorm(n)
  fit.ebmr = ebmr(X,y,tol=1e-10,admm=TRUE)
  fit.em = ridge_em1(y,X,1,1,200)
  expect_equal(fit.ebmr$residual_variance,fit.em$s2, tol=1e-3)
  expect_equal(fit.ebmr$g,fit.em$sb2/fit.em$s2, tol=1e-3)
  expect_true(all(diff(fit.ebmr$elbo)>=-1e8)) #check non-decreasing with tolerance
})
