#' @title Solve EBNV problem
#'
#' @description Functions for solving the Empirical Bayes Normal Variances (EBNV) problem for various
#' prior families
#'
#' @details These functions fit the model
#'
#' b_j~N(0,s2 w_j)
#'
#' w_j ~ g()
#'
#' by empirical
#' Bayes. They estimate the prior g, and (the inverse of the) posterior mean for each 1/w_j.
#'
#' @param b vector of n observations
#'
#' @param s2 a positive real number
#'
#' @param g (used only in some cases) an object with the same structure as the prior to be fit
#' used to specify things like the grid used for discrete priors.
#'
#' @return A list with elements
#'
#' \item{g}{A list containing the details of the estimated prior}
#'
#' \item{wbar}{A vector containing inverses of posterior mean for 1/w}
#'
#' @describeIn ebnv.np Solve EBNV problem with non-parametric prior
#' @export
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

#' @describeIn ebnv.np Solve EBNV problem with exponential prior
#' @inheritParams ebnv.np
#' @export
ebnv.exp = function(b,s2,g=NULL){
  w = 2*mean(abs(b))^2/s2
  wbar =  (1/sqrt(s2)) * (abs(b)*sqrt(w/2))
  g$w=w
  return(list(g=g,wbar=wbar))
}

#' @describeIn ebnv.np Solve EBNV problem with point mass prior
#' @inheritParams ebnv.np
#' @export
ebnv.pm = function(b,s2,g=NULL){
  w = mean(b^2)/s2
  wbar =  rep(w,length(b))
  g$w =w
  return(list(g=g,wbar=wbar))
}

#' @describeIn ebnv.np Solve EBNV problem with mixture of exponentials prior
#' @inheritParams ebnv.np
#' @export
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

