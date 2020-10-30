# Note that allows for non-zero prior mean -- ridge regression is
# usually 0 prior mean.
ridge = function(y, A, prior_variance, prior_mean = rep(0,ncol(A)),
                 residual_variance = 1){
  n = length(y)
  p = ncol(A)
  L = chol(crossprod(A) + (residual_variance/prior_variance)*diag(p))
  b = backsolve(L,crossprod(A,y) +
                (residual_variance/prior_variance) * prior_mean,
                transpose = TRUE)
  b = backsolve(L, b)
  #b = chol2inv(L) %*% (t(A) %*% y + (residual_variance/prior_variance)*prior_mean)
  return(drop(b))
}

prox_regression = function(x, t, y, A, residual_variance = 1){
  ridge(y,A,prior_variance = t,prior_mean = x,residual_variance)
}

# This is version of ridge above that takes SVD of A instead of A
# this effectively allows (A'A + kI)^{-1} to be computed efficiently for any k
ridge_svd = function(y, A.svd, prior_variance,
                     prior_mean = rep(0,nrow(A.svd$v)),
                     residual_variance = 1){
  n = length(y)
  p = nrow(A.svd$v)
  #bnew = (A.svd$v %*% (diag(A.svd$d) %*% (t(A.svd$u) %*% y)) + (residual_variance/prior_variance)*prior_mean)
  #bnew = t(A.svd$v) %*% bnew
  bnew = A.svd$d * crossprod(A.svd$u,y)
  bnew = bnew + crossprod(A.svd$v, prior_mean) *
                  (residual_variance/prior_variance)
  bnew = ((A.svd$d^2 + (residual_variance/prior_variance))^(-1)) * bnew
  bnew = A.svd$v %*% bnew
  return(drop(bnew))
}

prox_regression_svd = function(x, t, y, A.svd, residual_variance = 1){
  ridge_svd(y,A.svd,prior_variance = t,prior_mean = x,residual_variance)
}

# I use lamba = 1/2w where w is a vector of prior variances
prox_l2_het = function(x, t, lambda){
  prior_prec = 2*lambda # prior precision
  data_prec = 1/t
  return(x * data_prec/(data_prec + prior_prec))
}

update.mu.admm = function(fit){
  fit$mu = with(fit,prox_regression_svd(z - u,1/rho,y,X.svd))
  fit$z = with(fit,prox_l2_het(mu + u,1/rho,0.5/wbar))
  fit$u = with(fit,u + mu - z)
  return(fit)
}

