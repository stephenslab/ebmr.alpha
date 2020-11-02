# update wbar and g using ebnv (lasso) where g is exponential
ebmr.update.ebnv.lasso = function(fit){
  warning("ebnv lasso not tested; also need to work out ELBO computations")
  #need to work out the details...
  bfit = with(fit, sqrt(mu^2 + residual_variance * Sigma_diag)) #
  ebnv.res = ebnv_lasso(bfit, fit$sb2 * fit$residual_variance)

  fit$g$w = ebnv.res$w
  fit$wbar = ebnv.res$wbar

  return(fit)
}

# solve ebnv problem with exponential prior
# model is b_j \sim N(0, s2 w_j)
# w_j \sim g()
# where g is Exp(g$w), b and are known
# returns the mle \hat{g} and wbar which is inverse of posterior mean of 1/w
ebnv_lasso = function(b,s2){
  w = 2*mean(abs(b))^2/s2
  wbar =  (1/sqrt(s2)) * (abs(b)*sqrt(w/2))
  return(list(w=w,wbar=wbar))
}
