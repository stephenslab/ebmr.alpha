#' @description Updates parameters by fitting full GRR model in case n<p
#' @details Operates by performing a full SVD on Xtilde = XW^0.5, after which
#' other computations (EM algorithm to maximize sb2; s2, and posterior computations)
#' are cheap
#'
#' @param fit a ebmr fit
#'
#' @param tol a tolerance for the EM algorithm
#'
#' @param maxiter maximum number of iterations for EM algorithm
#'
#' @param compute_sigma_full set to true to compute full Sigma matrix for testing purposes
#'
#' @return an ebmr fit
ebmr.update.grr = function(fit, tol = 1e-3, maxiter = 1000, compute_Sigma_full=FALSE){

  if(fit$n>fit$p)
      warning("grr update is not designed for n>p; results may be unreliable")

  # svd computation
  Xtilde = t(fit$wbar^0.5 * t(fit$X))
  Xtilde.svd = svd(Xtilde)

  # maximum likelihood estimation
  ytilde = drop(t(Xtilde.svd$u) %*% fit$y)
  yt.em3 = ridge_indep_em3(ytilde, Xtilde.svd$d^2,tol, maxiter, s2.init = fit$residual_variance, sb2.init = fit$sb2)

  fit$residual_variance = yt.em3$s2
  fit$sb2 = yt.em3$sb2

  fit$Sigma_diag = Sigma1_diag_woodbury_svd(fit$wbar,Xtilde.svd,fit$sb2)

  if(compute_Sigma_full){
    fit$Sigma_full = Sigma1_woodbury_svd(fit$wbar,Xtilde.svd,fit$sb2)
  }

  fit$mu = mu1_woodbury_svd(fit$y, fit$wbar, Xtilde.svd, fit$sb2)

  fit$h2_term = -0.5*fit$p + 0.5*sum(log(fit$sb2*fit$wbar)) - 0.5*sum(log(1+fit$sb2*Xtilde.svd$d^2))

  return(fit)
}



# This code is based on https://stephens999.github.io/misc/ridge_em_svd.html
# It obtains maximum likelihood estimates for (s2,sb2)
# in the model y_i | theta_i \sim N(theta_i, s2), theta_i \sim N(0,sb2 s2 d2_i)
# where y,d2 are given and (s2,sb2) are to be estimated.
#
# However, internally it uses a different parameterization to make
# convergence faster (see the url above).
# Specifically it fits y_i | theta_i \sim N(sb theta_i,s^2), theta_i \sim N(0,l2 d2_i)
# so it returns s2hat = s2
# and sb2hat = sb2 *l2 / s2
ridge_indep_em3 = function(y, d2, tol=1e-3, maxiter=1000, s2.init = 1, sb2.init=1, l2.init=1){
  k = length(y)
  s2 = s2.init
  sb2 = sb2.init
  l2 = l2.init
  loglik = sum(dnorm(y,mean=0,sd = sqrt(sb2*l2*d2 + s2),log=TRUE))
  for(i in 1:maxiter){

    prior_var = d2*l2 # prior variance for theta
    data_var = s2/sb2 # variance of y/sb, which has mean theta
    post_var = 1/((1/prior_var) + (1/data_var)) #posterior variance of theta
    post_mean =  post_var * (1/data_var) * (y/sqrt(sb2)) # posterior mean of theta

    sb2 = (sum(y*post_mean)/sum(post_mean^2 + post_var))^2
    l2 = mean((post_mean^2 + post_var)/d2)

    r = y - sqrt(sb2) * post_mean # residuals
    s2 = mean(r^2 + sb2 * post_var)
    loglik= c(loglik,sum(dnorm(y,mean=0,sd = sqrt(sb2*l2*d2 + s2),log=TRUE)))
    if((loglik[i+1]-loglik[i])<tol) break
  }
  return(list(s2=s2,sb2=sb2*l2/s2,loglik=loglik))
}




