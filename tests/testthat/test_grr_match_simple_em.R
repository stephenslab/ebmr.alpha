test_that("grr results match simple EM",{
  set.seed(100)
  sd = 10
  n = 20
  p = 10
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  fit.ebmr = ebmr.init(X,y,1)
  fit.grr = ebmr.update.grr(fit.ebmr, compute_Sigma_full = TRUE, tol=1-15, maxiter = 10000)
  fit.em = ridge_em1(y,X,1,1,10000)

  expect_equal(compute.ridge.loglik(fit.grr), fit.em$loglik[length(fit.em$loglik)])
})
