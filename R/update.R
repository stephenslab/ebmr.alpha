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

# Update both mu and Sigma using direct approach with matrix inversion.
ebmr.update.mu.and.Sigma.direct = function(fit){
  fit = ebmr.update.Sigma.direct(fit)
  fit$mu = with(fit,drop(Sigma_full %*% Xty))
  return(fit)
}

# Return log-det of original matrix from cholesky decomposition.
chol2logdet = function(L){
  sum(log(diag(L)^2))
}

ebmr.update.Sigma.direct = function(fit){
  XtX = svd2XtX(fit$X.svd)
  LL = chol(XtX + diag(1/(fit$sb2*fit$wbar)))
  fit$Sigma_full = chol2inv(LL)
  fit$Sigma_diag = diag(fit$Sigma_full)
  fit$h2_term = -0.5*fit$p - 0.5*chol2logdet(LL) # log-determinant from cholesky
  return(fit)
}

# Update Sigma using the woodbury formula.
# Uses the approximation format X'X approximately L'L + D,
# where L is the first k right-singular vectors of X.
# If k = NULL, uses all of the svd computed at initialization,
# which is "exact" if this is the full svd.
ebmr.update.Sigma.woodbury = function(fit, k = NULL, compute_Sigma_full = FALSE){

  if(is.null(k)){
    k = length(fit$X.svd$d) # Use all the elements of svd.
  }

  # Compute the L and D from the (truncated) svd of X.
  L = with(fit$X.svd,d[1:k] * t(v[,1:k]))
  D = fit$d - colSums(L^2) # This should be zero for full SVD... but for future expansion

  ww = 1/(1/(fit$sb2*fit$wbar) + D) # effective weights
  Ltilde = t(t(L) * sqrt(ww)) # scale columns of L by sqrt(ww)

  H = diag(nrow(Ltilde)) + tcrossprod(Ltilde)
  H.chol = chol(H)
  H.inv = chol2inv(H.chol)

  # only necessary for testing
  if(compute_Sigma_full)
    fit$Sigma_full = diag(sqrt(ww)) %*%
      (diag(fit$p) - t(Ltilde) %*% H.inv %*% Ltilde) %*% diag(sqrt(ww))

  fit$Sigma_diag = ww * (1 - colSums(Ltilde * H.inv %*% Ltilde))

  # log-determinant from Cholesky.
  fit$h2_term = -0.5*fit$p + 0.5*sum(log(ww)) - 0.5*chol2logdet(H.chol)

  return(fit)
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
