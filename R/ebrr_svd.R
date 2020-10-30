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

  fit = update.mu.Sigma.rv.gscale(fit, tol, maxiter) # update all of mu, Sigma, residual variance and scaling of w

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

update.grr = function(fit, tol, maxiter){
  Xtilde = t(fit$wbar^0.5 * t(fit$X))
  Xtilde.svd = svd(Xtilde)
  ytilde = drop(t(fit$Xtilde.svd$u) %*% fit$y)
  yt.em3 = ridge_indep_em3(ytilde, Xtilde.svd$d^2,tol, maxiter, s2.init = fit$residual_variance, sb2.init = fit$sb2)
  fit$residual_variance = yt.em3$s2
  fit$sb2 = yt.em3$sb2

  fit$Sigma_diag = sigma_diag_from_w(Xtilde.svd,fit$wbar,fit$sb2)
  fit$mu1 = wbar * ytilde
}


