#' @description Fits Empirical Bayes multiple regression model
#' @details This is a new approach based on generalized ridge regression
#' @param X an n times p numeric matrix of covariates
#' @param y an n vector of responses
#' @param tol small real number controlling convergence tolerance; algorithm stops when elbo changes less than tol
#' @param maxiter integer indicating maximum number of iterations
#' @return an object of class "ebmr" that contains fit details
ebmr = function(X, y, tol=1e-10, admm = TRUE, maxiter = 1000){
  fit = ebmr.init(X,y)
  for(i in 1:maxiter){
    if(admm){
      fit = update.mu.admm(fit)
      fit = update.Sigma.woodbury(fit)
    } else {
      fit = update.mu.and.Sigma.direct(fit)
    }

    fit = update.residual_variance(fit)
    fit = update.w.and.g.ridge(fit)
    fit = update.elbo(fit)
    if(ebmr_get_elbodiff(fit) < tol) break
  }
  return(fit)
}

