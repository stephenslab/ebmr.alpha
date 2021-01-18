test_that("ebnv results make sense",{
  set.seed(1)
  n = 1e5
  eta = 2
  w = rexp(n,rate=eta)
  s2 = 9
  b = rnorm(n,0,sd= sqrt(s2*w))
  temp=ebnv.exp(b,s2)

  expect_equal(temp$g$w,1/eta,tol=1e-2)
  expect_equal(as.numeric(lm(w ~ temp$wbar)$coef[2]),1,tol=1e-2)

  temp=ebnv.np(b,s2,g=list(w=seq(1e-5,2,length=20)))

  expect_equal(sum(temp$g$mixprop * temp$g$w),1/eta,tol=1e-2)
  expect_equal(as.numeric(lm(w ~ temp$wbar)$coef[2]),1,tol=1e-1)

  # same but different eta and s2 just in case
  n = 1e5
  eta = 2
  w = rexp(n,rate=eta)
  s2 = 9
  b = rnorm(n,0,sd= sqrt(s2*w))
  temp=ebnv.exp(b,s2)

  expect_equal(temp$g$w,1/eta,tol=1e-2)
  expect_equal(as.numeric(lm(w ~ temp$wbar)$coef[2]),1,tol=1e-2)

  temp=ebnv.np(b,s2,g=list(w=seq(1e-5,2,length=20)))

  expect_equal(sum(temp$g$mixprop * temp$g$w),1/eta,tol=1e-2)
  expect_equal(as.numeric(lm(w ~ temp$wbar)$coef[2]),1,tol=1e-1)

})

test_that("ebnv.exp and ebnv.mix_exp logliks match when the mix solution is trivial",{
  set.seed(1)
  b = rexp(100)
  s2= 1
  res.exp = ebnv.exp(b,s2)
  w = c(0.01,res.exp$g$w, 1000) # set up grid so that only the middle one will be plausible
  res.exp_mix = ebnv.exp_mix(b,s2,list(w=w))
  expect_equal(res.exp_mix$loglik,res.exp$loglik, tol=1e-2)
})
