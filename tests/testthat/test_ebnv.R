test_that("ebnv results make sense",{
  set.seed(1)
  n = 1e5
  eta = 2
  w = rexp(n,rate=eta)
  s2 = 9
  b = rnorm(n,0,sd= sqrt(s2*w))
  temp=ebnv_lasso(b,s2)

  expect_equal(temp$w,1/eta,tol=1e-2)
  expect_equal(as.numeric(lm(w ~ temp$wbar)$coef[2]),1,tol=1e-2)

  temp=ebnv_np(b,s2,seq(1e-5,2,length=20))

  expect_equal(sum(temp$mixprop * temp$wgrid),1/eta,tol=1e-2)
  expect_equal(as.numeric(lm(w ~ temp$wbar)$coef[2]),1,tol=1e-1)

  # same but different eta and s2 just in case
  n = 1e5
  eta = 0.1
  w = rexp(n,rate=eta)
  s2 = 0.1
  b = rnorm(n,0,sd= sqrt(s2*w))
  temp=ebnv_lasso(b,s2)

  expect_equal(temp$w,1/eta,tol=1e-2)
  expect_equal(as.numeric(lm(w ~ temp$wbar)$coef[2]),1,tol=1e-2)

  temp=ebnv_np(b,s2,seq(1e-5,10/eta,length=20))

  expect_equal(sum(temp$mixprop * temp$wgrid),1/eta,tol=1e-1)
  expect_equal(as.numeric(lm(w ~ temp$wbar)$coef[2]),1,tol=1e-1)
})
