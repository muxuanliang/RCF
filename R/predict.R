#predict.rsfit get the predicted treatment
predict.rsfit <- function(fit, newx, lossType = 'logistic', type = 'transformed'){
  splitNumber <- length(fit)
  res_all <- sapply(fit, function(t){
    z = newx %*% t$fit$beta
    pre <- cbind(1,splines2::bSpline(z, knots=t$fit$knots, intercept = FALSE, Boundary.knots = t$fit$boundaryPoint))%*%t$fit$xi
    if(sum(abs(t$fit$xi))<1e-5){
      pre <- z
    }
    pre
  })
  if (type=='eta'){
    res_all <- link(res_all, lossType = lossType)
  }
  apply(res_all,1,mean)
}

# predict.rsfit gets the predicted treatment for a splitted fit
predict.rsfitSplit <- function(fit, newx, derivative = FALSE,lossType = 'logistic', type = 'value'){
  t <- fit
  z <- newx %*% t$beta
  spline <- splines2::bSpline(z, knots=t$knots, intercept = FALSE, Boundary.knots = t$boundaryPoint)
  pre <- cbind(1,spline)%*%t$xi
  if(sum(abs(t$xi))<1e-5){
    pre <- z
  }
  if(derivative){
    div <- nsd(z, knots = t$knots, intercept = FALSE, Boundary.knots = t$boundaryPoint)
    pre <- cbind(0, div)%*%t$xi
    if(sum(abs(t$xi))<1e-5){
      pre <- 1
    }
  }
  if (type=='eta'){
    pre <- link(pre, lossType = lossType)
  }
  pre
}
