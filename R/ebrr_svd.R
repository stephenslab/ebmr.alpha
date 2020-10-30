#' @title Empirical Bayes ridge regression by SVD
#' @description Fits Empirical Bayes ridge regression model by svd approach
#' @details This is based on my investigations at https://stephens999.github.io/misc/ridge_em_svd.html
#' It first performs an svd on X and then estimates the hyper-parameters by iterative EM
#' @param X an n times p numeric matrix of covariates
#' @param y an n vector of responses
#' @param tol small real number controlling convergence tolerance; algorithm stops when elbo changes less than tol
#' @param maxiter integer indicating maximum number of iterations
#' @return an object of class "ebmr" that contains fit details
#' @examples
#' set.seed(100)
#' n= 100
#' p = n
#' X = matrix(0,nrow=n,ncol=n)
#' for(i in 1:n){
#' X[i:n,i] = 1:(n-i+1)
#' }
#' btrue = rep(0,n)
#' btrue[40] = 8
#' btrue[41] = -8
#' y = X %*% btrue + rnorm(n)
#' plot(y,main="true (black); fitted values from RR (red)")
#' lines(X %*% btrue)
#' lines(X %*% coef(y.ebrr), col=2)
#'
#' @export
ebrr_svd = function(X, y, tol=1e-10, maxiter = 1000){
  fit = ebmr.init(X,y)

  fit = ebmr.update.grr(fit, tol, maxiter) # update all of mu, Sigma, residual variance and scaling of w

  return(fit)
}

