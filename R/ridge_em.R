# a simple implementation of an em algorithm for ridge regression
# y = Xb + e
# e \sim N(0,s2)
# b \sim N(0,sb2)
# uses EM to estimate s2 and sb2 as well as posterior mean E(b|y,X,s2,sb2)
ridge_em1 = function(y,X, s2,sb2, niter=10){
  XtX = t(X) %*% X
  Xty = t(X) %*% y
  yty = t(y) %*% y
  n = length(y)
  p = ncol(X)
  loglik = rep(0,niter)
  for(i in 1:niter){
    V = chol2inv(chol(XtX+ diag(s2/sb2,p)))

    SigmaY = sb2 *(X %*% t(X)) + diag(s2,n)
    loglik[i] = mvtnorm::dmvnorm(as.vector(y),sigma = SigmaY,log=TRUE)

    Sigma1 = s2*V  # posterior variance of b
    mu1 = as.vector(V %*% Xty) # posterior mean of b

    s2 = as.vector((yty + sum(diag(XtX %*% (mu1 %*% t(mu1) + Sigma1)))- 2*sum(Xty*mu1))/n)
    sb2 = mean(mu1^2+diag(Sigma1))

  }
  return(list(s2=s2,sb2=sb2,loglik=loglik,postmean=mu1))
}
