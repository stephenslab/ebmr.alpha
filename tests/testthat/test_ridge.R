test_that("ridge regression results match simple em", {
  set.seed(100)
  sd = 10
  n = 20
  p = 10
  X = matrix(rnorm(n*p),ncol=p)
  btrue = rnorm(p)
  y = X %*% btrue + sd*rnorm(n)
  fit.embr = ebmr(X,y,200)
  fit.em = ridge_em1(y,X,1,1,200)
  expect_equal(fit.embr$residual_variance,fit.em$s2, tol=1e-3)
  expect_equal(fit.embr$g,fit.em$sb2/fit.em$s2, tol=1e-3)
  expect_true(all(diff(fit.embr$elbo)>=-1e8)) #check non-decreasing with tolerance
})
