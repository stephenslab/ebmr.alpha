#' @title Solve EBNV problem with non-parametric prior
#' @description Solves the EBNV problem for non-parametric prior (discrete grid)
#' @details fits the model b_j \sim N(0,s2 w_j) and w_j \sim g() by empirical
#' Bayes, where g is non-parametric distribution approximated by a discrete distribution
#' on a grid
#' @param b n-vector of observations
#' @param s2 a scalar variance parameter
#' @param wgrid a discrete grid of non-negative values (unlike ash, 0 not allowed)
#' @return a list with elements mixprop, wgrid and wbar (inverse of posterior mean for 1/w)
ebnv_np = function(b, s2, wgrid){
  if(!all(wgrid>0)){
    stop("elements of wgrid must be non-negative")
  }
  std_obs = t(outer(b,sqrt(s2*wgrid),"/")) # standardized observations
  loglik.matrix = t(log(1/sqrt(s2*wgrid)) + dnorm(std_obs,log=TRUE))
  loglik.max = apply(loglik.matrix, 1, max)
  lik.matrix = exp(loglik.matrix-loglik.max)

  mixprop = mixsqp::mixsqp(lik.matrix,control = list(verbose=FALSE))$x

  postprob = t(mixprop * t(lik.matrix)) # likelihood * prior
  postprob = postprob/rowSums(postprob) # normalize

  wbar = 1/colSums((1/wgrid) * t(postprob))

  return(list(mixprop=mixprop, wgrid=wgrid, wbar = wbar))
}

# n=1000
# pi = rep(0.1,10)
# wgrid = seq(0.1,10,length=10)
# w = 10*runif(n)
# s2 = 1
# b = rnorm(n, 0, sqrt(s2*w))
# fit = ebnv_np(b,s2,wgrid)
# plot(b,sqrt(fit$wbar))
#
# b = rnorm(n, 0, sqrt(s2*0.001))
# fit = ebnv_np(b,s2,wgrid)
# plot(b,fit$wbar)


