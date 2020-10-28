#' @description Fits Empirical Bayes multiple regression model
#' @details This is a new approach based on generalized ridge regression
#' @param X an n times p numeric matrix of covariates
#' @param y an n vector of responses
#' @param tol small real number controlling convergence tolerance; algorithm stops when elbo changes less than tol
#' @param maxiter integer indicating maximum number of iterations
#' @param k integer indicating the rank of approximation to use for Sigma update (default of NULL uses full rank)
#' @param thin integer indicating how often to update Sigma (eg thin=10 means every 10 iterations).
#' @return an object of class "ebmr" that contains fit details
ebmr = function(X, y, tol=1e-10, admm = TRUE, maxiter = 1000, k=NULL, thin = 1){
  fit = ebmr.init(X,y)
  if(!is.null(k)){
    if(!admm){
      stop("k option currently only valid when admm=TRUE")
    }
  }
  for(i in 1:maxiter){
    if(admm){
      fit = update.mu.admm(fit)
      if((i-1) %% thin == 0){ # update Sigma only every thin iterations
        fit = update.Sigma.woodbury(fit, k=k)
      }
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

