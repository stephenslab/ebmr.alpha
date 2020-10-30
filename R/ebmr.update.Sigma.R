ebmr.update.Sigma = function(fit, method = c("woodbury","direct"), k = NULL, compute_Sigma_full = FALSE){
  method = match.arg(method)
  if(method == "direct")
    fit = ebmr.update.Sigma.direct(fit)
  else if(method == "woodbury"){
    fit = ebmr.update.Sigma.woodbury(fit, k, compute_Sigma_full)
  } else {
    stop("invalid method")
  }
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

  fit$h2_term = -0.5*fit$p + 0.5*sum(log(ww)) - 0.5*chol2logdet(H.chol)

  return(fit)
}
