ebmr.init = function(X,y){
  fit = list()
  fit$p = ncol(X)
  fit$n = nrow(X)
  fit$X = X
  fit$y = y
  fit$Xty = t(X) %*% y
  fit$XtX = t(X) %*% X
  fit$mu = rep(0,fit$p)
  fit$Sigma = matrix(0,nrow=fit$p,ncol=fit$p)
  fit$residual_variance = sd(y)^2
  fit$wbar = rep(1,fit$p)
  fit$g = 1
  fit$elbo = NULL
  fit$KLw = 0 # the KL from qW to prior g(W)
  return(fit)
}
