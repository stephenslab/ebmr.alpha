#' @description Fits Empirical Bayes ridge regression model by svd approach
#' @details This is based on my investigations at https://stephens999.github.io/misc/ridge_em_svd.html
#' It first performs an svd on X and then estimates the hyper-parameters by iterative EM
#' @param X an n times p numeric matrix of covariates
#' @param y an n vector of responses
#' @param tol small real number controlling convergence tolerance; algorithm stops when elbo changes less than tol
#' @param maxiter integer indicating maximum number of iterations
#' @return an object of class "ebmr" that contains fit details
ebrr_svd = function(X, y, tol=1e-10, maxiter = 1000){
  fit = ebmr.init(X,y)

  fit = update.grr.svd(fit, tol, maxiter) # update all of mu, Sigma, residual variance and scaling of w

  return(fit)
}

