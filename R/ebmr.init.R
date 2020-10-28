# fit is a list with elements
# Data elements:
# XtX
# Xty
# X
#
# Parameter elements:
# mu
# Sigma (diagonal elements of Sigma)
# wbar
# residual_variance
# g
#
# questions: should fit contain X or just fitted values Xmu?
# or maybe fit should contain svd for X? X=udv

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
  fit$KLw = 0 # the KL from qW to prior g(W)

  fit$h2_term = h2_func(fit)
  fit$elbo = elbo(fit)

  # for admm updates
  fit$z = fit$mu
  fit$u = fit$mu
  fit$rho = 1 # may want to change this

  class(fit) = "ebmr"
  return(fit)
}

ebmr.init.admm = function(f){

}
