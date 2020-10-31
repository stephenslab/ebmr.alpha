# update wbar and g using EBNV (non-parametric)
ebmr.update.ebnv.np = function(fit){
  if(length(fit$g$w)==1){
    stop("No grid set for w: add grid with ebmr.set.wgrid")
  }
  ebnv.res = ebnv_np(b = drop(sqrt(fit$mu^2 + fit$Sigma_diag)),
                     s2 = fit$sb2*fit$residual_variance, wgrid = fit$g$w)

  fit$g$mixprop = ebnv.res$mixprop
  fit$wbar = ebnv.res$wbar
  warning("need to implement KL and ElogW terms for this update")
  return(fit)
}
