# fit is a list with elements
# Data elements:
# Xty
# X.svd the svd of X
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
  fit$y = drop(y)
  fit$Xty = drop(crossprod(X,y))
  fit$X.svd = svd(X)
  fit$d = colSums(X^2)

  fit$mu = rep(0,fit$p)
  fit$Sigma_full = matrix(0,nrow = fit$p,ncol = fit$p) # for storing the full matrix; for testing only
  fit$Sigma_diag = rep(0,fit$p) # for storing the diagonal of Sigma

  fit$residual_variance = sd(y)^2
  fit$wbar = rep(1,fit$p)
  fit$g = 1
  fit$KLw = 0 # the KL from qW to prior g(W)
  fit$Elogw = 0
  
  fit$h2_term = -Inf # because Sigma is initialized at 0
  fit$elbo = elbo(fit)

  # For admm updates.
  fit$z = fit$mu
  fit$u = fit$mu
  fit$rho = 1 # may want to change this

  class(fit) = "ebmr"
  return(fit)
}


