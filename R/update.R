fitted.values = function(fit){
  fit$X %*% fit$mu
}

#' @description update residual variance by maximizing ELBO given all other parameters
update.residual_variance = function(fit){
  fit$residual_variance = (1/fit$n)*(sum((fit$y-fitted.values(fit))^2) + sum(fit$mu^2/fit$wbar))
  return(fit)
}

#' @description update both mu and Sigma using direct approach with matrix inversion
update.mu.and.Sigma.direct = function(fit){
  fit = update.Sigma.direct(fit)
  fit$mu = fit$Sigma_full %*% fit$Xty
  return(fit)
}

# return log-det of original matrix from cholesky decomposition
chol2logdet = function(L){
  sum(log(diag(L)^2))
}


update.Sigma.direct = function(fit){
  LL = chol(fit$XtX + diag(1/fit$wbar))
  fit$Sigma_full = chol2inv(LL)
  fit$Sigma_diag = diag(fit$Sigma_full)

  fit$h2_term =  -0.5*fit$p - 0.5 * chol2logdet(LL) # log-determinant from cholesky
  return(fit)
}

# update Sigma using the woodbury formula
# uses the approximation format X'X approximately L'L+D
# where L is the first k right-singular vectors of X
# if k=NULL uses all of the svd computed at initialization
# which is "exact" if this is the full svd
update.Sigma.woodbury = function(fit,k=NULL, compute_Sigma_full = FALSE){

  if(is.null(k)){k = length(fit$X.svd$d)} # use all the elements of svd

  # Compute the L and D from the (truncated) svd of X
  L = fit$X.svd$d[1:k] * t(fit$X.svd$v[,1:k])
  D = fit$d - colSums(L^2) # this should be 0 for full SVD... but for future expansion

  ww = 1/(1/fit$wbar + D) # effective weights
  Ltilde = t(t(L) * ww^0.5) # scale columns of L by ww^0.5

  H = diag(nrow(Ltilde)) + Ltilde %*% t(Ltilde)
  H.chol = chol(H)
  H.inv = chol2inv(H.chol)

  # only necessary for testing
  if(compute_Sigma_full)
    fit$Sigma_full = diag(ww^0.5) %*% (diag(fit$p) - t(Ltilde) %*% H.inv %*% Ltilde) %*% diag(ww^0.5)

  fit$Sigma_diag = ww * (1- colSums(Ltilde * H.inv %*% Ltilde))

  fit$h2_term =  -0.5*fit$p + 0.5* sum(log(ww)) - 0.5 * chol2logdet(H.chol) # log-determinant from cholesky
  return(fit)
}


#' @description update w and g in case where g is a point mass ("ridge regression)
update.w.and.g.ridge = function(fit){
  fit$g = mean(fit$mu^2 + fit$residual_variance * fit$Sigma_diag)/fit$residual_variance
  wbar_new = rep(fit$g,fit$p)
  fit = update.wbar(fit,wbar_new)
  fit$Elogw = rep(log(fit$g),fit$p) # the expected log term for elbo computation
  fit$KLw = 0 # the prior and posterior are the same for ridge

  return(fit)
}

# this function deals with the fact you need to update h2_term when updating wbar
# because the 0.5 tr(W^{-1} Sigma) term changes
update.wbar = function(fit,wbar_new){
  fit$h2_term = fit$h2_term + 0.5 * sum( ((1/fit$wbar)-(1/wbar_new)) * fit$Sigma_diag)
  fit$wbar = wbar_new
  return(fit)
}

update.elbo = function(fit){
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


# Sigma_full = chol2inv(chol(fit$XtX + diag(1/fit$wbar)))


# diag(Sigma_full) - Sigma_diag(fit$wbar,fit$X)

