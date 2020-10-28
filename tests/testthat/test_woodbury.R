test_that("woodbury update matches direct update", {
  set.seed(100)
  sd = 10
  n = 20
  p = 10
  X = matrix(rnorm(n*p),ncol=p)
  btrue = rnorm(p)
  y = X %*% btrue + sd*rnorm(n)
  fit.init = ebmr.init(X,y)
  fit1 = update.Sigma.direct(fit.init)
  fit2 = update.Sigma.woodbury(fit.init)
  expect_equal(fit1$Sigma,fit2$Sigma,tol=1e-8)
  expect_equal(fit1$h2_term,fit2$h2_term,tol=1e-8)
})
