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

  temp=ebnv.np(b,s2,g.init=list(w=seq(1e-5,2,length=20)))

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

  temp=ebnv.np(b,s2,g.init=list(w=seq(1e-5,2,length=20)))

  expect_equal(sum(temp$g$mixprop * temp$g$w),1/eta,tol=1e-2)
  expect_equal(as.numeric(lm(w ~ temp$wbar)$coef[2]),1,tol=1e-1)

})

test_that("ebnv.exp and ebnv.mix_exp logliks match when the mix solution is trivial",{
  set.seed(1)
  b = rexp(100)
  s2= 1
  res.exp = ebnv.exp(b,s2)
  w = c(0.01,res.exp$g$w, 1000) # set up grid so that only the middle one will be plausible
  res.exp_mix = ebnv.exp_mix(b,s2,g.init=list(w=w), update.mixprop= "mixsqp", update.w = "none")
  expect_equal(res.exp_mix$loglik,res.exp$loglik, tol=1e-2)
})


test_that("ebnv.mix_exp em updates are increasing log-likelihood",{
  set.seed(1)
  n = 1000
  w = c(rexp(n,rate = 1),rexp(n,rate=10))
  b = rnorm(2*n,0,sd=sqrt(w))

  g = list(mixprop=c(0.5,0.5), w=c(2,3))
  niter = 100
  loglik = rep(0,niter)
  for(i in 1:niter){
    res = ebnv.exp_mix(b,1,g,update.mixprop = "em",update.w="em")
    loglik[i] = res$loglik
    g = res$g
  }
  expect_true(all(diff(loglik)>=0))
})
