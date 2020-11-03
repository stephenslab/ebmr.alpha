# translate an exponential prior to an approximate non-parametric (discrete) equivalent
# currently just uses equal weights on quantiles
exp2np = function(g,grid.length=20){
  g$w = qexp(seq(0.01,0.99,length=grid.length),rate=1/g$w)
  g$mixprop = rep(1,grid.length)/grid.length
  return(g)
}
