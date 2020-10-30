# a simple implementation of an em algorithm for ridge regression
# y = Xb + e
# e \sim N(0,s2)
# b \sim N(0,sb2)
# uses EM to estimate s2 and sb2 as well as posterior mean E(b|y,X,s2,sb2)
#
#' @importFrom mvtnorm dmvnorm
ridge_em1 = function(y, X, s2, sb2, niter = 10){
  XtX = crossprod(X)
  Xty = drop(crossprod(X,y))
  yty = drop(crossprod(y))
  n = length(y)
  p = ncol(X)
  loglik = rep(0,niter)
  for(i in 1:niter){
    V = chol2inv(chol(XtX + diag(s2/sb2,p)))
    SigmaY = sb2*tcrossprod(X) + diag(s2,n)
    loglik[i] = dmvnorm(drop(y),sigma = SigmaY,log = TRUE)

    Sigma1 = s2*V  # posterior variance of b
    mu1 = drop(V %*% Xty) # posterior mean of b

    s2 = drop((yty + sum(diag(XtX %*% (tcrossprod(mu1) + Sigma1)))
                   - 2*sum(Xty*mu1))/n)
    sb2 = mean(mu1^2 + diag(Sigma1))
  }
  return(list(s2 = s2,sb2 = sb2,loglik = loglik,postmean = mu1))
}
