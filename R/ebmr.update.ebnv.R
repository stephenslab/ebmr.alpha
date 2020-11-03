# update wbar and g using ebnv (lasso) where g is exponential
ebmr.update.ebnv = function(fit,ebnv_fn){
  warning("need to work out ELBO computations for EBNV")

  bfit = with(fit, sqrt(mu^2 + residual_variance * Sigma_diag)) #
  ebnv.res = do.call(ebnv_fn,
                     list(b=bfit, s2=fit$sb2 * fit$residual_variance,
                          g = fit$g))

  fit$g = ebnv.res$g

  fit = ebmr.set.wbar(fit,ebnv.res$wbar)

  return(fit)
}


