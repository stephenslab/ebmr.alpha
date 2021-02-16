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
#' @param g.init (present only in some cases) an object with the same structure as the prior to be fit
#' used to initialize the algorithm (and so, implicitly, to sometimes specify parameters like the grid to be used)
#'
#' @return A list with elements
#'
#' \item{g}{A list containing the details of the estimated prior}
#'
#' \item{wbar}{A vector containing inverses of posterior mean for 1/w}
#'
#' \item{loglik}{The log-likelihood p(b | s2, ghat)}
#'
#'
#' @describeIn ebnv.pm Solve EBNV problem with point mass prior
#' @export
ebnv.pm = function(b,s2){
  w = mean(b^2)/s2
  wbar =  rep(w,length(b))
  g = list(mixprop = c(1), w = w)
  loglik = sum(dnorm(b,0,sqrt(s2*w), log=TRUE))
  return(list(g=g,wbar=wbar, loglik=loglik))
}


#' @describeIn ebnv.pm Solve EBNV problem with exponential prior
#' @inheritParams ebnv.pm
#' @export
ebnv.exp = function(b,s2){
  w = 2*mean(abs(b))^2/s2
  wbar =  (1/sqrt(s2)) * (abs(b)*sqrt(w/2))
  g = list(mixprop = c(1), w=w)

  # logllikelihood under double exponential on b
  loglik = sum(log(0.5) + dexp(b, rate= sqrt(2/(w*s2)), log=TRUE ))

  return(list(g=g,wbar=wbar, loglik = loglik))
}

compute_postprob = function(mixprop,lik.matrix){
  postprob = t(mixprop * t(lik.matrix)) # likelihood * prior
  postprob = postprob/rowSums(postprob) # normalize
  return(postprob)
}

#' @describeIn ebnv.pm Solve EBNV problem with non-parametric prior
#' @inheritParams ebnv.pm
#' @export
ebnv.np = function(b, s2, g.init, update.mixprop = c("em","mixsqp","none"), update.w = c("em","none")){
  g = g.init
  if(!all(g$w>0)){
    stop("elements of w grid must be positive")
  }
  update.mixprop = match.arg(update.mixprop)
  update.w = match.arg(update.w)


  std_obs = t(outer(b,sqrt(s2*g$w),"/")) # standardized observations
  loglik.matrix = t(log(1/sqrt(s2*g$w)) + dnorm(std_obs,log=TRUE))
  loglik.max = apply(loglik.matrix, 1, max)
  lik.matrix = exp(loglik.matrix-loglik.max)

  if(update.mixprop=="mixsqp"){
    res.mixsqp = mixsqp::mixsqp(lik.matrix,control = list(verbose=FALSE))
    g$mixprop = res.mixsqp$x
  } else if(update.mixprop=="em"){
    postprob = compute_postprob(g$mixprop,lik.matrix) # likelihood * prior
    g$mixprop = colMeans(postprob)
  }

  if(update.w == "em"){
    postprob = compute_postprob(g$mixprop,lik.matrix) # likelihood * prior

    b2bar = colSums(postprob*(b^2))/colSums(postprob) #weighted mean of abs(b)
    g$w = b2bar/s2

    std_obs = t(outer(b,sqrt(s2*g$w),"/")) # standardized observations
    loglik.matrix = t(log(1/sqrt(s2*g$w)) + dnorm(std_obs,log=TRUE))
    loglik.max = apply(loglik.matrix, 1, max)
    lik.matrix = exp(loglik.matrix-loglik.max)
  }

  postprob = compute_postprob(g$mixprop,lik.matrix) # likelihood * prior

  wbar = 1/colSums((1/g$w) * t(postprob))
  loglik = sum(log(colSums(t(lik.matrix)*g$mixprop)) + loglik.max)

  return(list(g=g, wbar = wbar, loglik=loglik))
}

#' @describeIn ebnv.pm Solve EBNV problem with non-parametric prior with EM update of both grid and mixture proportions
#' @inheritParams ebnv.pm
#' @export
ebnv.np.em = function(b, s2, g.init){
  ebnv.np(b, s2, g.init, "em", "em")
}

#' @describeIn ebnv.pm Solve EBNV problem with non-parametric prior with fixed grid
#' @inheritParams ebnv.pm
#' @export
ebnv.np.fixgrid = function(b, s2, g.init){
  ebnv.np(b, s2, g.init, "mixsqp", "none")
}


#' @describeIn ebnv.pm Solve EBNV problem with mixture of exponentials prior
#' @param update.mixprop string indicating how to estimate/update the mixture proportions; if "none" then mixture proportions are supplied by g$mixprop
#' @param update.w string indicating how to update w parameters; ; if "none" then mixture proportions are supplied by g$w
#' @export
ebnv.exp_mix = function(b, s2, g.init, update.mixprop = c("em","mixsqp","none"), update.w = c("em","none")){
  g=g.init
  if(!all(g$w>0)){
    stop("elements of w grid must be positive")
  }
  update.mixprop = match.arg(update.mixprop)
  update.w = match.arg(update.w)

  # first compute the likelihood matrix and posterior probabilities
  lambda = sqrt(2/(s2*g$w)) # the K rate parameters of implied double-exponential prior on b, one per gridpoint
  loglik.matrix = t(log(0.5) + log(lambda) - outer(lambda, abs(b), "*")) #n by K matrix
  loglik.max = apply(loglik.matrix, 1, max)
  lik.matrix = exp(loglik.matrix-loglik.max)

  if(update.mixprop=="mixsqp"){
    res.mixsqp = mixsqp::mixsqp(lik.matrix,control = list(verbose=FALSE))
    g$mixprop = res.mixsqp$x
  } else if(update.mixprop=="em"){
    postprob = compute_postprob(g$mixprop,lik.matrix) # likelihood * prior
    g$mixprop = colMeans(postprob)
  }

  if(update.w == "em"){
    postprob = compute_postprob(g$mixprop,lik.matrix) # likelihood * prior

    bbar = colSums(postprob*abs(b))/colSums(postprob) #weighted mean of abs(b)
    g$w = (2/s2) * bbar^2

    lambda = sqrt(2/(s2*g$w)) # equal to bbar, the K rate parameters of implied double-exponential prior on b, one per gridpoint
    loglik.matrix = t(log(0.5) + log(lambda) - outer(lambda, abs(b), "*")) #n by K matrix
    loglik.max = apply(loglik.matrix, 1, max)
    lik.matrix = exp(loglik.matrix-loglik.max)
  }

  postprob = compute_postprob(g$mixprop,lik.matrix) # likelihood * prior

  wbar = rowSums(outer(s2*abs(b)^(-1),lambda, "*") * postprob)^(-1)
  loglik = sum(log(colSums(t(lik.matrix)*g$mixprop)) + loglik.max)

  return(list(g=g, wbar = wbar, loglik = loglik))
}

#' @describeIn ebnv.pm Solve EBNV problem with mixture of exponentials prior with  EM update of both grid and mixture proportions
#' @inheritParams ebnv.pm
#' @export
ebnv.exp_mix.em = function(b, s2, g.init){
  ebnv.exp_mix(b, s2, g.init, "em", "em")
}

#' @describeIn ebnv.pm Solve EBNV problem with mixture of exponentials prior with fixed grid
#' @inheritParams ebnv.pm
#' @export
ebnv.exp_mix.fixgrid = function(b, s2, g.init){
  ebnv.exp_mix(b, s2, g.init, "mixsqp", "none")
}
