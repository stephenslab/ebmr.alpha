c_func = function(fit){
 -0.5*fit$n*log(2*pi*fit$residual_variance) - 0.5*sum(fit$Elogw)
}


h1_func = function(fit){
  -(0.5/fit$residual_variance) * (sum((fit$y - fitted.values(fit))^2) + sum(fit$mu^2/fit$wbar))
}

h2_func = function(fit){
  -0.5*sum( (fit$XtX + diag(1/fit$wbar)) * fit$Sigma ) + 0.5* as.numeric(determinant(fit$Sigma,logarithm=TRUE)$modulus)
}

elbo = function(fit){
  return(c_func(fit) + h1_func(fit) + h2_func(fit) + fit$KLw)
}
