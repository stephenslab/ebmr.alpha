ebmr.update.mu.direct = function(fit){
  if(is.null(fit$Sigma_full)){
    stop("update mu direct only works if Sigma_full is defined")
  }
  fit$mu = with(fit,drop(Sigma_full %*% Xty))
  return(fit)
}

fitted.values = function(fit){
  with(fit,drop(X %*% mu))
}

# Update residual variance by maximizing ELBO given all other
# parameters.
ebmr.update.residual_variance = function(fit){
  fit$residual_variance = with(fit,
    (1/n)*(sum((y - fitted.values(fit))^2) + sum(mu^2/(sb2*wbar))))
  return(fit)
}


# Return log-det of original matrix from cholesky decomposition.
chol2logdet = function(L){
  sum(log(diag(L)^2))
}


# Update w and g in case where g is a point mass ("ridge regression").
ebmr.update.ebnv.ridge = function(fit){
  fit$g = with(fit,mean(mu^2 + residual_variance*Sigma_diag)/residual_variance)
  wbar_new = rep(fit$g,fit$p)
  fit = update.wbar(fit,wbar_new)
  fit$Elogw = rep(log(fit$g),fit$p) # The expected log-term for elbo.
  fit$KLw = 0 # The prior and posterior are the same for ridge.
  return(fit)
}

# This function deals with the fact you need to update h2_term when
# updating wbar because the 0.5 tr(W^{-1} Sigma) term changes
update.wbar = function(fit, wbar_new){
  fit$h2_term = with(fit,
    h2_term + 0.5* (1/sb2) * sum(((1/wbar) - (1/wbar_new))*Sigma_diag))
  fit$wbar = wbar_new
  return(fit)
}

ebmr.update.elbo = function(fit){
  fit$elbo = c(fit$elbo,elbo(fit))
  return(fit)
}

# Sigma_diag = function(w, X){
#   Xtilde = X %*% diag(w^0.5)
#   n = nrow(X)
#   H = diag(n) + Xtilde %*% t(Xtilde)
#   Hinv = chol2inv(chol(H))
#   w * (1- colSums(Xtilde * Hinv %*% Xtilde))
# }
