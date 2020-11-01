# This code comes from my original implementation of veb lasso
#  https://stephens999.github.io/misc/lasso_em.html
# used as a check
calc_s2hat = function(y,X,XtX,EB,EB2){
  n = length(y)
  Xb = X %*% EB
  s2hat = as.numeric((1/n)* (t(y) %*% y - 2*sum(y*Xb) + sum(XtX * EB2)))
}

# this version if for fixed eta and s2
# if compute_mode=TRUE we have the regular LASSO
blasso_em = function(y,X,b.init,s2=1,eta=1,niter=100,compute_mode=FALSE){
  n = nrow(X)
  p = ncol(X)
  b = b.init # initiolize
  XtX = t(X) %*% X
  Xty = t(X) %*% y
  obj = rep(0,niter)
  EB = rep(1,p)
  varB = rep(1,p)

  for(i in 1:niter){
    W = as.vector(1/sqrt(varB + EB^2) * sqrt(2/eta))

    V = chol2inv(chol(XtX+ diag(s2*W)))
    Sigma1 = s2*V  # posterior variance of b
    varB = diag(Sigma1)
    if(compute_mode){
      varB = rep(0,p)
    }
    mu1 = as.vector(V %*% Xty) # posterior mean of b
    EB = mu1
  }
  return(list(bhat=EB))
}


# s2 is residual variance
# if compute_mode=TRUE we have the regular LASSO
blasso_veb = function(y,X,b.init,s2=1,eta=1,niter=100,update.eta=TRUE){
  n = nrow(X)
  p = ncol(X)
  b = b.init # initiolize
  XtX = t(X) %*% X
  Xty = t(X) %*% y
  obj = rep(0,niter)
  EB = rep(1,p)
  EB2 = diag(p)

  for(i in 1:niter){
    W = as.vector(1/sqrt(diag(EB2) + EB^2) * sqrt(2/eta))

    V = chol2inv(chol(XtX+ diag(s2*W)))
    Sigma1 = s2*V  # posterior variance of b
    varB = diag(Sigma1)

    mu1 = as.vector(V %*% Xty) # posterior mean of b
    EB = mu1
    EB2 = Sigma1 + outer(mu1,mu1)

    eta = mean(sqrt(diag(EB2))*sqrt(eta/2) + eta/2)
    s2 = calc_s2hat(y,X,XtX,EB,EB2)
  }
  return(list(bhat=EB,eta=eta,s2 = s2, W=W))
}
