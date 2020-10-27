ebmr = function(X, y, niter = 100){
  fit = ebmr.init(X,y)
  for(i in 1:niter){
    fit = update.mu.and.Sigma.full(fit)
    fit = update.residual_variance(fit)
    fit = update.w.and.g.ridge(fit)
    fit = update.elbo(fit)
  }
  return(fit)
}
