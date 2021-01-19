#' @title Fit Empirical Bayes Multiple Regression Model
#'
#' @description Fits multiple regression model using a new empirical Bayes variational
#' approach based on generalized ridge regression (GRR).
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
#' @param ebnv_fn Function for solving the Empirical Bayes Normal Variances problem.
#'
#' @param compute_mode Boolean indicating whether to use the mode (rather than the mean). If used
#' with ebnv_fn = ebnv.exp then this produces the Lasso solution rather than Bayesian Lasso.
#'
#' @return An object of class "ebmr" that contains fit details.
#'
#'
#' @examples
#' set.seed(100)
#' n = 50
#' p = 200
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' btrue = rep(0,p)
#' btrue[1] = 1
#' btrue[2] = -1
#' y = X %*% btrue + rnorm(n)
#'
#' y.fit.ebr = ebmr(X,y, maxiter = 200, ebnv_fn = ebnv.pm)
#' y.fit.eblasso = ebmr(X,y, maxiter = 200, ebnv_fn = ebnv.exp)
#' y.fit.eblasso.mode = ebmr(X,y, maxiter = 200, ebnv_fn = ebnv.exp, compute_mode=TRUE)
#'
#' y.fit.ebash = ebmr.alpha:::ebmr.set.prior(y.fit.eblasso,ebmr.alpha:::exp2np(y.fit.eblasso$g))
#' y.fit.ebash = ebmr.update(y.fit.ebash, maxiter = 200, ebnv_fn = ebnv.np)
#'
#' plot(btrue, coef(y.fit.ebr))
#' plot(btrue, coef(y.fit.eblasso), col=2)
#' plot(btrue, coef(y.fit.eblasso.mode), col=3)
#' plot(btrue, coef(y.fit.ebash), col=4)
#' @export
ebmr = function (X, y, tol = 1e-8, maxiter = 1000, ebnv_fn = ebnv.exp, compute_mode = FALSE){
  fit = ebmr.init(X,y)
  fit = ebmr.update(fit, tol, maxiter, ebnv_fn, compute_mode)
  return(fit)
}

#' @describeIn ebmr Updates a previous EBMR fit, usually using a new prior family
#' @param fit The previous EBMR fit
#' @export
ebmr.update = function (fit, tol = 1e-8, maxiter = 1000, ebnv_fn = ebnv.exp, compute_mode = FALSE, update_residual_variance = TRUE){

  for(i in 1:maxiter){

    fit = ebmr.update.grr.svd(fit, compute_Sigma_diag = !compute_mode, update_residual_variance = update_residual_variance)
    #fit = ebmr.scale.sb2(fit)

    fit = ebmr.update.ebnv(fit,ebnv_fn)

    fit = ebmr.update.elbo(fit)
    if(abs(ebmr_get_elbodiff(fit)) < tol) # use absolute until the elbo computation is correct
      break
  }
  return(fit)
}


# @param admm If \code{admm = TRUE}, use the ADMM updates for
#   updating the regression coefficients.
# @param k Integer indicating the rank of the approximation to use
#   for Sigma update (default \code{NULL} uses full rank)
#
# @param thin Integer indicating how often to update Sigma (e.g.,
#   \code{thin = 10} means update every 10 iterations).
#
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

