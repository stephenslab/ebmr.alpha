
# updates mu (finds posterior mean of b) in a GRR
# using ADMM
# convergence tolerance is on the difference between mu and z
ebmr.update.mu.admm = function(fit, maxiter = 1000, tol=1e-3){
  for(i in 1:maxiter){
    fit$mu = with(fit,prox_regression_svd(z - u,1/rho,y,X.svd))
    fit$z = with(fit,prox_l2_het(mu + u,1/rho,0.5/(sb2*wbar)))
    fit$u = with(fit,u + mu - z)
    if(all(abs(fit$mu-fit$z)<tol))
      break
  }
  return(fit)
}

# updates both mu and Sigma by direct method (may as well update Sigma as well
# as mu because the matrix inversion is the hard part)
ebmr.update.mu.Sigma.direct = function(fit){
  fit = ebmr.update.Sigma.direct(fit)
  fit$mu = with(fit,drop(Sigma_full %*% Xty))
  return(fit)
}
