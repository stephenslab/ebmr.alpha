c_func = function(fit){
 -0.5*fit$n*log(2*pi*fit$residual_variance) - 0.5*sum(fit$Elogw)
}


h1_func = function(fit){
  -(0.5/fit$residual_variance) * (sum((fit$y - fitted.values(fit))^2) + sum(fit$mu^2/fit$wbar))
}

# this function is not directly used in the algorithm - only for testing
h2_func = function(fit){
  XtX = svd2XtX(X.svd)
  -0.5*sum( (XtX + diag(1/fit$wbar)) * fit$Sigma_full ) + 0.5* as.numeric(determinant(fit$Sigma_full,logarithm=TRUE)$modulus)
}

elbo = function(fit){
  return(c_func(fit) + h1_func(fit) + fit$h2_term + fit$KLw)
}

# compute X'X from the svd of X
svd2XtX = function(X.svd){
  return(X.svd$v %*% (X.svd$d^2 * t(X.svd$v)))
}
