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
#
#' @title Initialize ebmr fit
#' @param X a numeric n times p matrix
#' @param y an n vector
#' @param sb2 scalar parameter value
#' @importFrom stats sd
#' @export
ebmr.init = function(X,y,sb2=1){
  fit = list()
  fit$p = ncol(X)
  fit$n = nrow(X)
  fit$X = X
  fit$y = drop(y)
  fit$Xty = drop(crossprod(X,y))
  fit$X.svd = svd(X)
  fit$d = colSums(X^2)

  fit$mu = rep(0,fit$p)
  fit$Sigma_full = NULL # for storing the full matrix; for testing only
  fit$Sigma_diag = rep(0,fit$p) # for storing the diagonal of Sigma

  fit$residual_variance = sd(y)^2
  fit$wbar = rep(1,fit$p)
  fit$sb2 = sb2

  # initialize to simple ridge case
  fit$g = list(mixprop = c(1), w= c(1))
  #fit$KLw = 0 # the KL from qW to prior g(W)
  #fit$Elogw = rep(log(fit$g$w),fit$p) # The expected log-term for elbo.

  fit$logdet_KL_term = -0.5*sum(log(fit$g$w))

  fit$h2_term = -Inf # because Sigma is initialized at 0
  fit$elbo = elbo(fit)

  # For admm updates.
  fit$z = fit$mu
  fit$u = fit$mu
  fit$rho = 1 # may want to change this

  class(fit) = "ebmr"
  return(fit)
}


