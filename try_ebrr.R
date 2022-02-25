set.seed(100)
n = 100
p = n
X = matrix(0,nrow=n,ncol=n)
for(i in 1:n){
  X[i:n,i] = 1:(n-i+1)
}
btrue = rep(0,n)
btrue[40] = 8
btrue[41] = -8
s2 = .4
y = X %*% btrue + s2*rnorm(n)

library(ebmr.alpha)

y.fit.ebr = ebmr(X,y, maxiter = 200, ebnv_fn = ebnv.pm)
y.fit.eblasso.fixs2 = ebmr.update(y.fit.ebr, maxiter = 200, ebnv_fn = ebnv.exp, update_residual_variance = FALSE)
y.fit.eblasso = ebmr.update(y.fit.eblasso.fixs2, maxiter = 200, ebnv_fn = ebnv.exp, update_residual_variance = TRUE)

y.fit.eblasso.mode = ebmr.update(y.fit.eblasso, maxiter = 200, ebnv_fn = ebnv.exp, compute_mode=TRUE)

plot(y,main="true (black)")
lines(X %*% btrue)
lines(X %*% coef(y.fit.ebr), col=2)
lines(X %*% coef(y.fit.eblasso.fixs2), col=3)

lines(X %*% coef(y.fit.eblasso),col=4)
lines(X %*% coef(y.fit.eblasso.mode),col=5)

ebnv.exp_mix.fix = function(b,s2){
  K=20 # grid size
  gridval = 10*y.fit.eblasso$g$w*(2^(0.05*((1:K)-1)) - 1)^2 + 1e-5
  return(ebnv.exp_mix(b,s2,g=list(mixprop=rep(1/K,K), w= gridval), update.mixprop="none", update.w="none"))
}

ebnv.exp_mix.mixsqp = function(b,s2,g){
  ebnv.exp_mix(b,s2,g,update.mixprop = "mixsqp", update.w="none")
}

ebnv.exp_mix.em = function(b,s2,g){
  ebnv.exp_mix(b,s2,g,update.mixprop = "em", update.w = "em")
}

y.fit.eb.expmix.fix = ebmr.update(y.fit.eblasso,ebnv_fn = ebnv.exp_mix.fix)
y.fit.eb.expmix = ebmr.update(y.fit.eb.expmix.fix, maxiter = 200, ebnv_fn = ebnv.exp_mix.mixsqp)

y.fit.eb.expmix.em = ebmr.update(y.fit.eb.expmix.fix, maxiter = 200, ebnv_fn = ebnv.exp_mix.em)



#y.fit.ebash = ebmr.alpha:::ebmr.set.prior(y.fit.eblasso.mode,ebmr.alpha:::exp2np(y.fit.eblasso$g))
#y.fit.ebash = ebmr.update(y.fit.ebash, maxiter = 200, ebnv_fn = ebnv.np)


plot(y,main="true (black)")
lines(X %*% btrue)
lines(X %*% coef(y.fit.eb.expmix.fix),col=2)
lines(X %*% coef(y.fit.eb.expmix),col=3)
lines(X %*% coef(y.fit.eb.expmix.em),col=4)


lines(X %*% coef(y.fit.ebash),col=7)


set.seed(100)
n = 100
p = 100
X = matrix(rnorm(n*p),nrow=n,ncol=p)

btrue = rep(0,n)
btrue[40] = 8
btrue[41] = -8
y = X %*% btrue + rnorm(n)

y.fit.ebr = ebmr(X,y, maxiter = 200, ebnv_fn = ebnv.pm)
y.fit.eblasso = ebmr.update(y.fit.ebr, maxiter = 200, ebnv_fn = ebnv.exp)

y.fit.ebash = ebmr.alpha:::ebmr.set.prior(y.fit.eblasso,ebmr.alpha:::exp2np(y.fit.eblasso$g))
y.fit.ebash = ebmr.update(y.fit.ebash, maxiter = 200, ebnv_fn = ebnv.np)

plot(btrue, coef(y.fit.ebr))
plot(btrue, coef(y.fit.eblasso), col=2)
plot(btrue, coef(y.fit.ebash), col=3)



set.seed(100)
n = 50
p = 200
X = matrix(rnorm(n*p),nrow=n,ncol=p)

btrue = rep(0,p)
btrue[1] = 1
btrue[2] = -1
y = X %*% btrue + rnorm(n)

y.fit.ebr = ebmr(X,y, maxiter = 200, ebnv_fn = ebnv.pm)
y.fit.eblasso = ebmr(X,y, maxiter = 200, ebnv_fn = ebnv.exp)

plot(y.fit.ebr$mu)
plot(y.fit.eblasso$mu)

sqrt(mean((y.fit.ebr$mu- btrue)^2))
sqrt(mean((y.fit.eblasso$mu- btrue)^2))


y.fit.ebr$residual_variance
y.fit.eblasso$residual_variance

plot(y.fit.eblasso$elbo)
plot(y.fit.ebr$elbo,col=2)

my_ebnv = function(b,s2,g=NULL){
  ebnv.exp_mix(b,s2,g=list(w=seq(0.05,20,length=20)))
}
y.fit.expmix = ebmr(X,y, maxiter = 200, ebnv_fn = my_ebnv)




#null simulation
set.seed(1)
n = 50
p = 200
X = matrix(rnorm(n*p),nrow=n,ncol=p)

btrue = rep(0,p)
y = X %*% btrue + rnorm(n)

y.fit.ebr = ebmr(X,y, maxiter = 200, ebnv_fn = ebnv.pm)
y.fit.eblasso = ebmr.update(y.fit.ebr, maxiter = 200, ebnv_fn = ebnv.exp)

plot(y.fit.ebr$mu)
plot(y.fit.eblasso$mu)

sqrt(mean((y.fit.ebr$mu- btrue)^2))
sqrt(mean((y.fit.eblasso$mu- btrue)^2))




Xw = t(t(X) * y.fit.eblasso$wbar^0.5)
temp = ebmr(Xw,y)


res.mrash = mr.ash.alpha::mr.ash(X,y)

y.fit.ebash =
y.fit.ebash = ebmr.set.prior(y.fit.eblasso,ebmr.alpha:::exp2np(y.fit.eblasso$g))
y.fit.ebash = ebmr.update(y.fit.ebash, maxiter = 200, ebnv_fn = ebnv.np)

plot(btrue, coef(y.fit.ebr))
plot(btrue, coef(y.fit.eblasso), col=2)
plot(btrue, coef(y.fit.ebash), col=3)






y.blasso = ebmr.alpha:::blasso_em(y,X,drop(y.fit$mu),
                     s2 = y.fit$residual_variance,
                     eta = y.fit$g$w * y.fit$sb2 * y.fit$residual_variance,
                     niter = 500)

plot(y.blasso$bhat, y.fit$mu)

set.seed(1)
n=1000
p = 100
X = matrix(rnorm(n*p), nrow=n)
sb2 = 4
btrue = rnorm(p,sd = sqrt(sb2))
#btrue = rexp(p, rate=1)
s2 = 1
y = X %*% cbind(btrue) + rnorm(n,sd=sqrt(s2))


y.init = ebmr.init(X,y)
y.grr = ebmr.update.grr(y.init, tol=1e-10)



y.fit = ebmr(X,y,maxiter=200,tol=1e-3)
y.fit$residual_variance

y.blasso = blasso_veb(y,X,btrue,niter = 200)

y.blasso$s2

y.blasso$eta

y.fit$residual_variance * y.fit$g$w * y.fit$sb2

