compute.ridge.loglik = function(fit){
  SigmaY = with(fit, residual_variance * sb2 * g$w *(X %*% t(X)) + diag(residual_variance,n))

  loglik = mvtnorm::dmvnorm(as.vector(y),sigma = SigmaY,log=TRUE)
  return(loglik)
}
