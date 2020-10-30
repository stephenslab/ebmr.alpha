context("woodbury")

test_that("woodbury update matches direct update",{
  set.seed(100)
  sd = 10
  n = 20
  p = 10
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  fit.init = ebmr.init(X,y,sb2=2)
  fit1 = ebmr.update.Sigma.direct(fit.init)
  fit2 = ebmr.update.Sigma.woodbury(fit.init, compute_Sigma_full = TRUE)
  expect_equal(fit1$Sigma_full,fit2$Sigma_full,scale = 1,tol = 1e-8)
  expect_equal(fit1$Sigma_diag,fit2$Sigma_diag,scale = 1,tol = 1e-8)
  expect_equal(fit1$h2_term,fit2$h2_term,scale = 1,tol = 1e-8)
})
