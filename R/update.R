fitted.values = function(fit){
  fit$X %*% fit$mu
}

#' @description update residual variance by maximizing ELBO given all other parameters
update.residual_variance = function(fit){
  fit$residual_variance = (1/fit$n)*(sum((fit$y-fitted.values(fit))^2) + sum(fit$mu^2/fit$wbar))
  return(fit)
}

#' @description update both mu and Sigma using direct approach with matrix inversion
update.mu.and.Sigma.full = function(fit){
  fit = update.Sigma.full(fit)
  fit$mu = fit$Sigma %*% fit$Xty
  return(fit)
}

update.Sigma.full = function(fit){
  fit$Sigma = chol2inv(chol(fit$XtX + diag(1/fit$wbar)))
  return(fit)
}


#' @description update w and g in case where g is a point mass ("ridge regression)
update.w.and.g.ridge = function(fit){
  fit$g = mean(fit$mu^2 + fit$residual_variance * diag(fit$Sigma))/fit$residual_variance
  fit$wbar = rep(fit$g,fit$p)
  fit$Elogw = rep(log(fit$g),fit$p) # the expected log term for elbo computation
  fit$KLw = 0 # the prior and posterior are the same for ridge
  return(fit)
}

update.elbo = function(fit){
  fit$elbo = c(fit$elbo,elbo(fit))
  return(fit)
}
