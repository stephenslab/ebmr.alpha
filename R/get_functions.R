ebmr_get_elbodiff = function(fit){
  niter = length(fit$elbo)
  return(fit$elbo[niter] - fit$elbo[niter-1])
}
