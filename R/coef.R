#' @title extract regression coefficients from ebmr fit
#' @param object an ebmr fit
#' @return a p vector of estimated regression coefficients
#' @export
coef.ebmr = function(object, ...){
  object$mu
}



