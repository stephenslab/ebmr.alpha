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

# solve ebnv problem with exponential prior
# model is b_j \sim N(0, s2 w_j)
# w_j \sim g()
# where g is Exp(mean = g$w), b and are known
# returns the mle \hat{g} and wbar which is inverse of posterior mean of 1/w
ebnv_exp = function(b,s2){
  w = 2*mean(abs(b))^2/s2
  wbar =  (1/sqrt(s2)) * (abs(b)*sqrt(w/2))
  return(list(w=w,wbar=wbar))
}


# solve ebnv problem with point mass prior
# model is b_j \sim N(0, s2 w_j)
# w_j \sim g()
# where g is delta(g$w), b and s2 are known
# returns the mle \hat{g} and wbar which is inverse of posterior mean of 1/w
ebnv_pm = function(b,s2){
  w = mean(b^2)/s2
  wbar =  rep(w,length(b))
  return(list(w=w,wbar=wbar))
}

