# update wbar and g using ebnv (lasso) where g is exponential
ebmr.update.ebnv = function(fit,ebnv_fn){
  warning("need to work out ELBO computations for EBNV")

  bfit = with(fit, sqrt(mu^2 + residual_variance * Sigma_diag)) #
  s2 = with(fit, sb2 * residual_variance)

  ebnv.res = do.call(ebnv_fn,
                     list(b=bfit, s2=s2,
                          g = fit$g))

  fit$g = ebnv.res$g
  fit = ebmr.set.wbar(fit,ebnv.res$wbar)
  fit$logdet_KL_term = ebnv.res$loglik +
    (fit$p/2)*log(2*pi*s2) + (0.5/s2) * sum(bfit^2/fit$wbar)

  return(fit)
}


