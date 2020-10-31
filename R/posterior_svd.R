# compute posterior quantites for GRR using SVD of Xtilde = XW^0.5


# computes sigma1 by direct inversion, for comparison
Sigma1_direct = function(w,X,sb2=1){
  chol2inv(chol(diag(1/(w*sb2))+t(X) %*% X))
}

# compute sigma1 by Woodbury from SVD
Sigma1_woodbury_svd = function(w,Xt.svd,sb2=1){
  u = Xt.svd$u
  v = Xt.svd$v
  d = Xt.svd$d
  p = nrow(v)
  dt = sb2 * d^2/(1+sb2*d^2)
  sb2*diag(w^0.5) %*% (diag(1,p) - v %*%  diag(dt) %*% t(v)) %*% diag(w^0.5)
}

# compute diagonal elements of Sigma1 by Woodbury from SVD
Sigma1_diag_woodbury_svd = function(w,Xt.svd,sb2){
  u = Xt.svd$u
  v = Xt.svd$v
  d = Xt.svd$d
  dt = sb2 * d^2/(1+sb2*d^2)
  sb2 * w * (1 - colSums(t(v^2)*dt))
}


mu1_direct = function(y,w,X,sb2=1){
  drop(chol2inv(chol(diag(1/(w*sb2))+t(X) %*% X)) %*% t(X) %*% y)
}

mu1_woodbury_svd = function(y,w,Xt.svd,sb2=1){
  u = Xt.svd$u
  v = Xt.svd$v
  d = Xt.svd$d
  dt = sb2 * d^2/(1+sb2*d^2)
  #sb2 * diag(w^0.5) %*% v %*% diag(1-dt) %*% diag(d) %*% t(u) %*% y

  #diag(w^0.5) %*% v %*% diag(dt/d) %*% t(u) %*% y
  drop(w^0.5 * ( v %*% ((dt/d) * (t(u) %*% y))))
}

