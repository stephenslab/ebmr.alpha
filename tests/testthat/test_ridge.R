context("ridge")

test_that("grr results do not change after simple updates",{
  #grr should not change with further updates as they are optimal for the grr
  set.seed(100)
  sd = 10
  n = 20
  p = 20
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  fit.ebmr = ebmr.init(X,y,1)
  fit.grr = ebmr.update.grr.svd(fit.ebmr, compute_Sigma_full = TRUE, tol=1-15)
  fit.grr = ebmr.update.elbo(fit.grr)

  fit.rr = ebmr.update.Sigma.woodbury(fit.grr, compute_Sigma_full = TRUE) # should not move as it should be optimal

  expect_equal(fit.grr,fit.rr,tol=1e-8)

  suppressWarnings(fit.rr <- ebmr.update.ebnv(fit.rr, ebnv_fn = ebnv.pm))
  fit.grr$logdet_KL_term = fit.rr$logdet_KL_term # force these to be equal since fit.grr does not include that term

  expect_equal(fit.grr,fit.rr,tol=1e-8)

  fit.rr = ebmr.update.residual_variance(fit.rr)
  expect_equal(fit.grr,fit.rr,tol=1e-8)

  fit.rr = ebmr.update.mu.Sigma.direct(fit.rr)
  expect_equal(fit.grr,fit.rr,tol=1e-8)

  #scale by sb2 and then refit
  fit.grr = ebmr.scale.sb2(fit.grr)
  fit.rr = ebmr.update.Sigma.woodbury(fit.grr, compute_Sigma_full = TRUE) # should not move as it should be optimal
  expect_equal(fit.rr$Sigma_diag,fit.grr$Sigma_diag,tol=1e-8)

  fit.rr = ebmr.update.residual_variance(fit.rr)
  expect_equal(fit.grr$residual_variance,fit.rr$residual_variance,tol=1e-8)

  fit.rr = ebmr.update.mu.Sigma.direct(fit.rr)
  expect_equal(fit.grr$mu,fit.rr$mu,tol=1e-8)
})


test_that("grr results match simple em",{
  set.seed(100)
  sd = 10
  n = 20
  p = 20
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  fit = ebmr.init(X,y)
  fit.grr = ebmr.update.grr.svd(fit,tol = 1e-10,maxiter=1000)
  fit.em = ridge_em1(y,X,1,1,1000)
  expect_equal(fit.grr$residual_variance,fit.em$s2, tol=1e-3)
  expect_equal(fit.grr$sb2*fit.grr$g$w,fit.em$sb2/fit.em$s2, tol=1e-3)
  expect_true(all(diff(fit.grr$elbo)>=-1e8)) #check non-decreasing with tolerance
})

test_that("mu results from direct and admm updates match grr",{
  set.seed(100)
  sd = 10
  n = 20
  p = 20
  X = matrix(rnorm(n*p),ncol = p)
  btrue = rnorm(p)
  y = drop(X %*% btrue + sd*rnorm(n))
  fit = ebmr.init(X,y)
  fit.grr = ebmr.update.grr.svd(fit,tol = 1e-10,maxiter=1000)
  fit.admm = ebmr.update.mu.admm(fit.grr,tol = 1e-10,maxiter=1000)
  fit.direct = ebmr.update.mu.Sigma.direct(fit.grr)

  expect_equal(fit.admm$mu, fit.grr$mu, tol=1e-5)
  expect_equal(fit.direct$mu, fit.grr$mu, tol=1e-5)

  fit$sb2 = 1
  fit$g$w = fit.grr$g$w * fit.grr$sb2
  fit$residual_variance = fit.grr$residual_variance
  fit$wbar = rep(fit$g$w,p)

  fit.admm = ebmr.update.mu.admm(fit,tol = 1e-10,maxiter=1000)
  fit.direct = ebmr.update.mu.Sigma.direct(fit)
  expect_equal(fit.admm$mu, fit.grr$mu, tol=1e-5)
  expect_equal(fit.direct$mu, fit.grr$mu, tol=1e-5)

  # Test for non-constant wbar
  fit = ebmr.set.wbar(fit, runif(p))
  fit.grr = ebmr.update.grr.svd(fit,tol = 1e-10,maxiter=1000)
  fit.admm = ebmr.update.mu.admm(fit.grr,tol = 1e-10,maxiter=1000)
  fit.direct = ebmr.update.mu.Sigma.direct(fit.grr)

  expect_equal(fit.admm$mu, fit.grr$mu, tol=1e-5)
  expect_equal(fit.direct$mu, fit.grr$mu, tol=1e-5)


})
