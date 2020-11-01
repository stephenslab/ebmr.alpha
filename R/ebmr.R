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
#' @param admm If \code{admm = TRUE}, use the ADMM updates for
#'   updating the regression coefficients.
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
 ebmr = function (X, y, tol = 1e-10, maxiter = 1000){ # admm = FALSE, k=NULL){
  fit = ebmr.init(X,y)

  for(i in 1:maxiter){

    fit = ebmr.update.grr(fit)

    # if(admm){
    #   fit = ebmr.update.grr.admm(fit, k=k, maxiter = 1)
    # } else {
    #   fit = ebmr.update.grr.direct(fit)
    # }
    #fit = ebmr.update.ebnv.ridge(fit)
    fit = ebmr.update.ebnv.lasso(fit)

    fit = ebmr.update.elbo(fit)
    if(abs(ebmr_get_elbodiff(fit)) < tol) # use absolute until the elbo computation is correct
      break
  }
  return(fit)
}

# # fits grr model by admm updates on mu
# # updating Sigma every thin iterations
# # and using rank k approximation
# ebmr.update.grr.admm = function(fit, tol = 1e-10, maxiter = 1000, k = NULL, thin = 1){
#   warning("admm updates not fully tested; in particular tolerance check on ELBO may not be good idea")
#   for(i in 1:maxiter){
#     fit = ebmr.update.mu.admm(fit)
#     fit = ebmr.update.residual_variance(fit)
#     # Update Sigma only every "thin" iterations.
#     if((i-1) %% thin == 0){
#       fit = ebmr.update.Sigma.woodbury(fit,k = k)
#     }
#     fit = ebmr.update.elbo(fit)
#     if(abs(ebmr_get_elbodiff(fit)) < tol) # check absolute value as admm may not be monotonic?
#       break
#   }
#
#   return(fit)
# }


# fits grr model by direct (naive) updates on mu, Sigma
ebmr.update.grr.direct = function(fit, tol = 1e-10, maxiter = 1000){
  for(i in 1:maxiter){
    fit = ebmr.update.Sigma.direct(fit)
    fit = ebmr.update.mu.direct(fit)
    fit = ebmr.update.residual_variance(fit)
    fit = ebmr.update.elbo(fit)
    if(ebmr_get_elbodiff(fit) < tol)
      break
  }
  return(fit)
}

