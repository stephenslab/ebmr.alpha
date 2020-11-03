#' @export
ebmr.set.prior = function(fit,g){
  fit$g = g
  return(fit)
}

# This function deals with the fact you need to update h2_term when
# updating wbar because the 0.5 tr(W^{-1} Sigma) term changes
#' @export
ebmr.set.wbar = function(fit, wbar_new){
  fit$h2_term = with(fit,
                     h2_term + 0.5* (1/sb2) * sum(((1/wbar) - (1/wbar_new))*Sigma_diag))
  fit$wbar = wbar_new
  return(fit)
}
