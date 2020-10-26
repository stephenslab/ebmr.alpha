# fit is a list with elements
# Data elements:
# XtX
# Xty
# X
#
# Parameter elements:
# mu
# Sigma (diagonal elements of Sigma)
# w
# residual_variance
# g
#
# questions: should fit contain X or just fitted values Xmu?
# or maybe fit should contain svd for X? X=udv

fitted.values = function(fit){
  fit$X %*% fit$mu
}

h1 = function(fit){
  -(0.5/fit$residual_variance) * (sum((fit$y - fitted.values(y))^2) + sum(fit$mu^2/w))
}



elbo.ridge = function(fit){

}

#' @description update residual variance by maximizing ELBO given all other parameters
update.residual_variance = function(fit){
  n = nrow(fit$X)
  fit$residual_variance = (1/n)*(sum((fit$y-fitted.values(fit))^2) + sum(fit$mu^2/fit$w))
  return(fit)
}

#' @description update both mu and Sigma using direct approach with matrix inversion
update.mu.and.Sigma.full = function(fit){
  Sigma = chol2inv(chol(fit$XtX + diag(1/fit$w)))
  fit$mu = Sigma %*% fit$Xty
  fit$Sigma = diag(Sigma)
  return(fit)
}


#' @description update w and g in case where g is a point mass ("ridge regression)
update.w.and.g.ridge = function(fit){
  fit$g = mean(fit$mu^2 + fit$residual_variance * fit$Sigma)/fit$residual_variance
  fit$w = rep(fit$g,length(fit$Xty))
  return(fit)
}
