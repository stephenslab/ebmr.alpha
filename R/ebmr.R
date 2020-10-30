#' @title Fit Empirical Bayes Multiple Regression Model
#' 
#' @description This is a new approach based on generalized ridge
#'   regression.
#' 
#' @param X An n times p numeric matrix of covariates.
#' 
#' @param y An n vector of responses.
#' 
#' @param tol Small real number controlling convergence tolerance;
#'   algorithm stops when elbo changes less than \code{tol}.
#' 
#' @param maxiter Integer indicating maximum number of iterations.
#' 
#' @param k Integer indicating the rank of the approximation to use
#'   for Sigma update (default \code{NULL} uses full rank)
#' 
#' @param thin Integer indicating how often to update Sigma (e.g.,
#'   \code{thin = 10} means update every 10 iterations).
#' 
#' @return An object of class "ebmr" that contains fit details.
#'
#' @export
#' 
ebmr = function (X, y, tol = 1e-10, admm = TRUE, maxiter = 1000,
                 k = NULL, thin = 1){
  fit = ebmr.init(X,y)
  if(!is.null(k)){
    if(!admm){
      stop("k option currently only valid when admm = TRUE")
    }
  }
  for(i in 1:maxiter){
    if(admm){
      fit = update.mu.admm(fit)

      # Update Sigma only every "thin" iterations.
      if((i-1) %% thin == 0){ 
        fit = update.Sigma.woodbury(fit,k = k)
      }
    } else {
      fit = update.mu.and.Sigma.direct(fit)
    }

    fit = update.residual_variance(fit)
    fit = update.w.and.g.ridge(fit)
    fit = update.elbo(fit)
    if(ebmr_get_elbodiff(fit) < tol)
      break
  }
  return(fit)
}
