ebmr.init = function(X,y){
  p = ncol(X)
  fit = list()
  fit$X = X
  fit$y = y
  fit$Xty = t(X) %*% y
  fit$XtX = t(X) %*% X
  fit$mu = rep(0,p)
  fit$Sigma = rep(0,p)
  fit$residual_variance = sd(y)^2
  fit$w = rep(1,p)
  fit$g = 1
  return(fit)
}
