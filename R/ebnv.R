#' @title Solve EBNV problem with non-parametric prior
#' @description Solves the EBNV problem for non-parametric prior (discrete grid)
#' @details fits the model b_j \sim N(0,s2 w_j) and w_j \sim g() by empirical
#' Bayes, where g is non-parametric distribution approximated by a discrete distribution
#' on a grid
#' @param b n-vector of observations
#' @param s2 a scalar variance parameter
#' @param g a list with elements mixprop and w, the latter being a non-negative grid of possible values (unlike ash, 0 not allowed)
#' @return a list with elements g and wbar (inverse of posterior mean for 1/w)
ebnv.np = function(b, s2, g){
  if(!all(g$w>0)){
    stop("elements of wgrid must be non-negative")
  }
  std_obs = t(outer(b,sqrt(s2*g$w),"/")) # standardized observations
  loglik.matrix = t(log(1/sqrt(s2*g$w)) + dnorm(std_obs,log=TRUE))
  loglik.max = apply(loglik.matrix, 1, max)
  lik.matrix = exp(loglik.matrix-loglik.max)

  mixprop = mixsqp::mixsqp(lik.matrix,control = list(verbose=FALSE))$x

  postprob = t(mixprop * t(lik.matrix)) # likelihood * prior
  postprob = postprob/rowSums(postprob) # normalize

  wbar = 1/colSums((1/g$w) * t(postprob))
  g$mixprop = mixprop

  return(list(g=g, wbar = wbar))
}

# solve ebnv problem with exponential prior
# model is b_j \sim N(0, s2 w_j)
# w_j \sim g()
# where g is Exp(mean = g$w), b and are known
# returns the mle \hat{g} and wbar which is inverse of posterior mean of 1/w
ebnv.exp = function(b,s2,g=NULL){
  w = 2*mean(abs(b))^2/s2
  wbar =  (1/sqrt(s2)) * (abs(b)*sqrt(w/2))
  g$w=w
  return(list(g=g,wbar=wbar))
}


# solve ebnv problem with point mass prior
# model is b_j \sim N(0, s2 w_j)
# w_j \sim g()
# where g is delta(g$w), b and s2 are known
# returns the mle \hat{g} and wbar which is inverse of posterior mean of 1/w
ebnv.pm = function(b,s2,g=NULL){
  w = mean(b^2)/s2
  wbar =  rep(w,length(b))
  g$w =w
  return(list(g=g,wbar=wbar))
}

# mix ebnv for mixture of exponentials
# g is mixture of exponentials, with means given by grid

ebnv.exp_mix = function(b, s2, g){
  if(!all(g$w>0)){
    stop("elements of wgrid must be non-negative")
  }
  lambda = sqrt(2/(s2*g$w)) # the K rate parameters of double-exponential, one per gridpoint
  loglik.matrix = t(log(lambda) - outer(lambda, abs(b), "*")) #n by K matrix
  loglik.max = apply(loglik.matrix, 1, max)
  lik.matrix = exp(loglik.matrix-loglik.max)

  mixprop = mixsqp::mixsqp(lik.matrix,control = list(verbose=FALSE))$x

  postprob = t(mixprop * t(lik.matrix)) # likelihood * prior
  postprob = postprob/rowSums(postprob) # normalize

  wbar = rowSums(outer(s2*abs(b)^(-1),lambda, "*") * postprob)^(-1)

  g$mixprop = mixprop

  return(list(g=g, wbar = wbar))
}

