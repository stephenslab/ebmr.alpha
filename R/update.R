
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
