# Update w and g in case where g is a point mass ("ridge regression").
ebmr.update.ebnv.ridge = function(fit){
  fit$g = with(fit,mean(mu^2 + residual_variance*Sigma_diag)/(sb2*residual_variance))
  wbar_new = rep(fit$g,fit$p)
  fit = update.wbar(fit,wbar_new)
  fit$Elogw = rep(log(fit$g),fit$p) # The expected log-term for elbo.
  fit$KLw = 0 # The prior and posterior are the same for ridge.
  return(fit)
}
