# update wbar and g using EBNV (non-parametric)
update.ebnv.np = function(fit){
  ebnv.res = ebnv_np(b = sqrt(fit$mu^2 + fit$Sigma_diag),
                     s2 = fit$sb2*fit$residual_variance, wgrid = fit$g$w)

  fit$g$mixprop = ebnv.res$mixprop
  fit$wbar = ebnv.res$wbar

}
