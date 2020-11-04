logdet = function (x)
  as.numeric(determinant(x,logarithm = TRUE)$modulus)

c_func = function(fit){
  with(fit,p/2 -0.5*n*log(2*pi*residual_variance) - (p/2)*log(sb2) - 0.5*sum(Elogw))
}

h1_func = function(fit){
  with(fit,-(0.5/residual_variance) * (sum((y - fitted.values(fit))^2)
                                       + sum(mu^2/(sb2*wbar))))
}

# This function is not directly used in the algorithm - only for testing.
h2_func = function(fit){
  XtX = svd2XtX(fit$X.svd)
  return(with(fit,
         -0.5*sum((XtX + diag(1/(sb2*wbar))) * Sigma_full) + 0.5*logdet(Sigma_full)))
}

elbo = function(fit){
  return(c_func(fit) + h1_func(fit) + fit$h2_term + fit$KLw)
}

# Compute X'X from the svd of X.
svd2XtX = function(X.svd){
  return(with(X.svd,v %*% (d^2 * t(v))))
}
