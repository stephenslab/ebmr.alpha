# compute posterior quantites for GRR using SVD of Xtilde = XW^0.5


# computes sigma1 by direct inversion, for comparison
Sigma1_direct = function(w,X,sb2=1){
  chol2inv(chol(diag(1/(w*sb2))+t(X) %*% X))
}

# compute sigma1 by Woodbury from SVD
Sigma1_wood_svd = function(w,Xt.svd,sb2=1){
  u = Xt.svd$u
  v = Xt.svd$v
  d = Xt.svd$d
  dt = sb2 * d^2/(1+sb2*d^2)
  sb2*diag(w^0.5) %*% (diag(1,p) - v %*%  diag(dt) %*% t(v)) %*% diag(w^0.5)
}

# compute diagonal elements of Sigma1 by Woodbury from SVD
Sigma1_diag_wood_svd = function(w,Xt.svd,sb2){
  u = Xt.svd$u
  v = Xt.svd$v
  d = Xt.svd$d
  dt = sb2 * d^2/(1+sb2*d^2)
  sb2 * w * (1 - colSums(t(v^2)*dt))
}

# testing
n = 10
p = 20
X = matrix(rnorm(n*p),nrow=n)
w = rnorm(p)^2
Xtilde = t(w^0.5 * t(X))
Xtilde.svd = svd(Xtilde)
sb2 = 0.5

all.equal(Sigma1_direct(w,X,sb2),Sigma1_wood_svd(w,Xtilde.svd,sb2))
all.equal(diag(Sigma1_direct(w,X,sb2)),Sigma1_diag_wood_svd(w,Xtilde.svd,sb2))

sb2 = 2

all.equal(Sigma1_direct(w,X,sb2),Sigma1_wood_svd(w,Xtilde.svd,sb2))
all.equal(diag(Sigma1_direct(w,X,sb2)),Sigma1_diag_wood_svd(w,Xtilde.svd,sb2))

mu1_direct = function(y,w,X,sb2=1){
  chol2inv(chol(diag(1/(w*sb2))+t(X) %*% X)) %*% t(X) %*% y
}

mu1_wood_svd = function(y,w,Xt.svd,sb2=1){
  u = Xt.svd$u
  v = Xt.svd$v
  d = Xt.svd$d
  dt = sb2 * d^2/(1+sb2*d^2)
  #sb2 * diag(w^0.5) %*% v %*% diag(1-dt) %*% diag(d) %*% t(u) %*% y

  diag(w^0.5) %*% v %*% diag(dt/d) %*% t(u) %*% y
}

y = rnorm(n)
mu1_direct(y,w,X,sb2) - mu1_wood_svd(y,w,Xt.svd,sb2)
