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
  fit$mu = fit$Sigma %*% fit$Xty
  return(fit)
}

# return log-det of original matrix from cholesky decomposition
chol2logdet = function(L){
  sum(log(diag(L)^2))
}


update.Sigma.direct = function(fit){
  LL = chol(fit$XtX + diag(1/fit$wbar))
  fit$Sigma = chol2inv(LL)
  fit$h2_term =  -0.5*fit$p - 0.5 * chol2logdet(LL) # log-determinant from cholesky
  return(fit)
}

# update Sigma using the woodbury formula
update.Sigma.woodbury = function(fit){
  Xtilde = fit$X %*% diag(fit$wbar^0.5)
  H = diag(fit$n) + Xtilde %*% t(Xtilde)
  H.chol = chol(H)

  H.inv = chol2inv(H.chol)
  fit$Sigma = diag(fit$wbar^0.5) %*% (diag(fit$p) - t(Xtilde) %*% H.inv %*% Xtilde) %*% diag(fit$wbar^0.5)

  fit$h2_term =  -0.5*fit$p + 0.5* sum(log(fit$wbar)) - 0.5 * chol2logdet(H.chol) # log-determinant from cholesky
  return(fit)
}


#' @description update w and g in case where g is a point mass ("ridge regression)
update.w.and.g.ridge = function(fit){
  fit$g = mean(fit$mu^2 + fit$residual_variance * diag(fit$Sigma))/fit$residual_variance
  wbar_new = rep(fit$g,fit$p)
  fit = update.wbar(fit,wbar_new)
  fit$Elogw = rep(log(fit$g),fit$p) # the expected log term for elbo computation
  fit$KLw = 0 # the prior and posterior are the same for ridge

  return(fit)
}

# this function deals with the fact you need to update h2_term when updating wbar
# because the 0.5 tr(W^{-1} Sigma) term changes
update.wbar = function(fit,wbar_new){
  fit$h2_term = fit$h2_term + 0.5 * sum( ((1/fit$wbar)-(1/wbar_new)) * diag(fit$Sigma))
  fit$wbar = wbar_new
  return(fit)
}

update.elbo = function(fit){
  fit$elbo = c(fit$elbo,elbo(fit))
  return(fit)
}

Sigma_diag = function(w, X){
  Xtilde = X %*% diag(w^0.5)
  n = nrow(X)
  H = diag(n) + Xtilde %*% t(Xtilde)
  Hinv = chol2inv(chol(H))
  w * (1- colSums(Xtilde * Hinv %*% Xtilde))
}


# Sigma_full = chol2inv(chol(fit$XtX + diag(1/fit$wbar)))


# diag(Sigma_full) - Sigma_diag(fit$wbar,fit$X)

