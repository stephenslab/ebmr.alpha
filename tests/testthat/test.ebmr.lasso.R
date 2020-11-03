test_that("eb lasso results match previous veb_lasso implementation",{
  set.seed(1)
  n=1000
  p = 100
  X = matrix(rnorm(n*p), nrow=n)
  sb2 = 4
  #btrue = rnorm(p,sd = sqrt(sb2))
  btrue = rexp(p, rate=1)
  s2 = 1
  y = X %*% cbind(btrue) + rnorm(n,sd=sqrt(s2))

  suppressWarnings(y.fit <- ebmr(X,y,maxiter=200,tol=1e-3))
  y.fit$residual_variance

  y.blasso = blasso_veb(y,X,btrue,niter = 200)

  expect_equal(y.blasso$s2,y.fit$residual_variance,tol=1e-3)
  expect_equal(y.blasso$eta, y.fit$residual_variance * y.fit$g$w * y.fit$sb2, tol=1e-3)
  expect_equal(y.blasso$bhat, y.fit$mu,tol=1e-3)
}
)
