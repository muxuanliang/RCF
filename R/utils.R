# utils
#
# This file codes each step in the algorithm proposed in the paper.
# In addition, it also provides defined link function can be used.
#

# ks gets the kernel estimation
ks <- function(xx, yy, xx.test){
  nobs <- nrow(as.matrix(xx))
  nvars <- ncol(as.matrix(xx))
  hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+4)))
  wm <- function(t){
    if (ncol(as.matrix(xx))==1){
      weight <- exp(-0.5 * (as.numeric(t)-xx)^2/(hopt^2)) * hopt
    } else {
      weight <- apply(xx,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(xx))})
    }
    weighted.mean(yy, weight)
  }
  if (nrow(as.matrix(xx.test))==1) {
    yy.test <- wm(xx.test)
  } else {
    if (ncol((as.matrix(xx.test)))==1){
      yy.test <- sapply(as.matrix(xx.test), function(t){
        wm(t)
      })
    } else {
      yy.test <- apply(as.matrix(xx.test),1,function(t){
        wm(t)
      })
    }
  }
  yy.test
}

# model.gam gets the gam regression function given data
model.gam <- function(data){
  p <- dim(data$predictor)[2]
  expr <- "mgcv::gam(outcome~"
  for (i in 1:(p-1)){
    expr <- paste0(expr, "s(predictor[,",i,"])+")
  }
  expr <- paste0(expr, "s(predictor[,",p,"]), data = data,method ='REML')")
  expr
}

# getOutcomeModel gets the outcome regression model
getOutcomeModel <- function(data, method=c('lm', 'glmnet', 'kernel', 'others'), splitIndex, Formula = NULL, predictAll = FALSE, screeningMethod="SIRS", outcomeScreeningFamily='Gaussian'){
  sampleSplitIndex <- splitIndex$nuisance
  p <- dim(data$predictor)[2]
  size <- dim(data$predictor)[1]
  fit <- NULL
  supp <- NULL
  dataPredict <- NULL
  dataPredict$control=data$predictor[sampleSplitIndex,]
  dataPredict$treatment=data$predictor[sampleSplitIndex,]
  if (predictAll){
    dataPredict$control=data$predictor
    dataPredict$treatment=data$predictor
  }
  prediction <- NULL
  dataControl <- list(predictor=data$predictor[(!sampleSplitIndex) & (data$treatment==FALSE),], outcome=data$outcome[(!sampleSplitIndex) & (data$treatment==FALSE)])
  dataTreatment <- list(predictor=data$predictor[(!sampleSplitIndex) & (data$treatment==TRUE),], outcome=data$outcome[(!sampleSplitIndex) & (data$treatment==TRUE)])
  if ((0.05*size >= p) | (p <= 5)){
    supp$control <- supp$treatment <- rep(TRUE, times = p)
    if ((method != 'lm')&&(method != 'glmnet')){
      ans1 <- screening(data$predictor[data$treatment==FALSE,], data$outcome[data$treatment==FALSE], method = screeningMethod, family = outcomeScreeningFamily)
      ans2 <- screening(data$predictor[data$treatment==TRUE,], data$outcome[data$treatment==TRUE], method = screeningMethod, family = outcomeScreeningFamily)
      if (screeningMethod == 'glmnet'){
        supp$control <- ans1
        supp$treatment <- ans2
      } else {
        supp$control <- supp$treatment <- (ans1 <= floor(p/2))|(ans2 <= floor(p/2))
      }
      dataControl$predictor <- dataControl$predictor[,supp$control]
      dataTreatment$predictor <- dataTreatment$predictor[,supp$treatment]
    }
  }
  if ((0.05*size < p) || (method == 'glmnet')) {
    fit$control <- glmnet::cv.glmnet(x = dataControl$predictor, y = dataControl$outcome)
    fit$treatment <- glmnet::cv.glmnet(x = dataTreatment$predictor, y = dataTreatment$outcome)
    supp$control <- abs(fit$control$glmnet.fit$beta[,fit$control$glmnet.fit$lambda==fit$control$lambda.min])>0
    supp$treatment <- abs(fit$treatment$glmnet.fit$beta[,fit$treatment$glmnet.fit$lambda==fit$treatment$lambda.min])>0
    if ((method != 'lm')&&(method != 'glmnet')){
      ans1 <- screening(data$predictor[data$treatment==FALSE,], data$outcome[data$treatment==FALSE], method = screeningMethod, family = outcomeScreeningFamily)
      ans2 <- screening(data$predictor[data$treatment==TRUE,], data$outcome[data$treatment==TRUE], method = screeningMethod, family = outcomeScreeningFamily)
      if (screeningMethod == 'glmnet'){
        supp$control <- ans1
        supp$treatment <- ans2
      } else {
        supp$control <- supp$treatment <- (ans1 <= 5)|(ans2 <= 5)
      }
    }
    dataControl$predictor <- dataControl$predictor[,supp$control]
    dataTreatment$predictor <- dataTreatment$predictor[,supp$treatment]
    prediction$control <- predict(fit$control, newx = dataPredict$control, s=fit$control$lambda.min)
    prediction$treatment <- predict(fit$treatment, newx = dataPredict$treatment, s=fit$treatment$lambda.min)
  }

  if (is.null(Formula)){
    Formula <- function(support){
      expr <- (outcome ~ predictor)
      if (sum(support)==0){
        expr <- (outcome ~ 1)
      }
      expr
    }
  }
  fit <- NULL
  dataPredict <- NULL
  dataPredict$control=data$predictor[sampleSplitIndex,supp$control]
  dataPredict$treatment=data$predictor[sampleSplitIndex,supp$treatment]
  if (predictAll){
    dataPredict$control=data$predictor[,supp$control]
    dataPredict$treatment=data$predictor[,supp$treatment]
  }
  if ((method == 'lm')||(method == 'glmnet')){
    if (sum(supp$control) > 0){
      fit$control <- lm(Formula(supp$control), data = dataControl)
      prediction$control <- predict(fit$control, newdata = list(predictor=dataPredict$control))
    }
    if (sum(supp$treatment) > 0){
      fit$treatment <- lm(Formula(supp$treatment), data = dataTreatment)
      prediction$treatment <- predict(fit$treatment, newdata = list(predictor=dataPredict$treatment))
    }
  } else if (method == 'kernel') {
    if (sum(supp$control) > 0){
      prediction$control <- ks(dataControl$predictor, dataControl$outcome, dataPredict$control)
    }
    if (sum(supp$treatment) > 0){
      prediction$treatment <- ks(dataTreatment$predictor, dataTreatment$outcome, dataPredict$treatment)
    }
  } else {
    if (sum(supp$control) > 0){
      fit$control <- eval(parse(text=model.gam(dataControl)))
      prediction$control <- predict(fit$control, newdata = list(predictor=dataPredict$control))
    }
    if (sum(supp$treatment) > 0){
      fit$treatment <- eval(parse(text=model.gam(dataTreatment)))
      prediction$treatment <- predict(fit$treatment, newdata = list(predictor=dataPredict$treatment))
    }
  }
  prediction
}

# getPropensityModel gets the estimation of the propensity model based on selected method
getPropensityModel <- function(data, method=c('lm', 'glmnet', 'kernel'), splitIndex, Formula = NULL, predictAll = FALSE, screeningMethod="SIRS"){
  sampleSplitIndex <- splitIndex$nuisance
  p <- dim(data$predictor)[2]
  size <- dim(data$predictor)[1]
  fit <- NULL
  supp <- NULL
  dataPredict <- NULL
  dataPredict=data$predictor[sampleSplitIndex,]
  if (predictAll){
    dataPredict=data$predictor
  }
  prediction <- NULL
  dataTrain <- list(predictor=data$predictor[(!sampleSplitIndex),], treatment=data$treatment[(!sampleSplitIndex)])
  if (0.05*size >= p){
    supp$control <- supp$treatment <- rep(TRUE, times = p)
    if ((method != 'lm')&&(method != 'glmnet')){
      ans <- screening(data$predictor, data$treatment, method = screeningMethod, family = 'binomial')
      if (screeningMethod == 'glmnet'){
        supp <- ans
      } else {
        supp <- (ans <= p/2)
      }
    }
    dataTrain$predictor = dataTrain$predictor[,supp]
  }
  if ((0.05*size < p) || (method == 'glmnet')) {
    fit <- glmnet::cv.glmnet(x = dataTrain$predictor, y = dataTrain$treatment, family='binomial')
    supp <- abs(fit$glmnet.fit$beta[,fit$glmnet.fit$lambda==fit$lambda.min])>0
    if ((method != 'lm')&&(method != 'glmnet')){
      ans <- screening(data$predictor, data$treatment, method = screeningMethod, family = 'binomial')
      if (screeningMethod == 'glmnet'){
        supp <- ans
      } else {
        supp <- (ans <= 5)
      }
    }
    dataTrain$predictor <- dataTrain$predictor[,supp]
    prediction <- predict(fit, newx = dataPredict, type='response', s=fit$lambda.min)
  }

  if (is.null(Formula)){
    Formula <- function(support){
      expr <- (treatment ~ predictor)
      if (sum(support)==0){
        expr <- (treatment ~ 1)
      }
      expr
    }
  }
  fit <- NULL
  dataPredict <- NULL
  dataPredict=data$predictor[sampleSplitIndex,supp]
  if (predictAll){
    dataPredict=data$predictor[,supp]
  }
  if ((method == 'lm')||(method == 'glmnet')){
    if (sum(supp) > 0){
      fit <- glm(Formula(supp), family=binomial, data = dataTrain)
      prediction <- predict(fit, newdata = list(predictor=dataPredict), type="response")
    }
  } else if (method == 'kernel') {
    if (sum(supp) > 0){
      prediction <- ks(dataTrain$predictor, dataTrain$treatment, dataPredict)
      prediction <- (prediction > 0.9) * 0.9 + (prediction < 0.1) * 0.1 + (prediction < 0.9) * (prediction > 0.1) * prediction
    }
  }
  prediction
}

# screening gets the screened variable for high dimeniosnal data
screening <- function(x, y, method='glmnet', family='Gaussian'){
  var <- apply(x, 2, sd)
  supp <- order(var, decreasing = TRUE)
  if (method=='glmnet'){
    fit <- glmnet::cv.glmnet(x, y, family = family)
    coef <- fit$glmnet.fit$beta[,fit$lambda==fit$lambda.min]
    supp <- (abs(coef)>0)
  } else {
    fit <- VariableScreening::screenIID(x, y, method=method)
    supp <- fit$rank
  }
  supp
}

# loss reutrn the value or the derivative of the chosen type of link
loss <- function(input, type = c('logistic', 'exponential'), order = 'original'){
  return(switch(order, 'original'=switch(type, 'logistic'=exp(-input)/2, 'exponential'=exp(-2*input)/4),
         'derivative'=switch(type, 'logistic'=-exp(-input)/2, 'exponential'=-exp(-2*input)/2),
         'hessian'=switch(type, 'logistic'=exp(-input)/2, 'exponential'=exp(-2*input))))
}

# link returns the eta from eta_trans
link <- function(eta_trans, lossType='logistic'){
  return(switch(lossType, 'logistic'=(exp(eta_trans/2)-exp(-eta_trans/2))/(exp(eta_trans/2)+exp(-eta_trans/2)),
                'exponential'=(exp(eta_trans)-exp(-eta_trans))/(exp(eta_trans)+exp(-eta_trans))))
}

# nsd get the derivative of the ns
nsd <- function(x, knots = NULL, intercept = FALSE, Boundary.knots = range(x)){
  div <- sapply(x,function(t){
    if((t+1e-6)>=Boundary.knots[2]){
      res <- (splines2::bSpline(t, knots=knots, intercept=intercept, Boundary.knots =Boundary.knots)-splines2::bSpline(t-1e-6, knots=knots, intercept=intercept, Boundary.knots =Boundary.knots))/(1e-6)
    } else if ((t-1e-6)<=Boundary.knots[1]){
      res <- (splines2::bSpline(t+1e-6, knots=knots, intercept=intercept, Boundary.knots =Boundary.knots)-splines2::bSpline(t, knots=knots, intercept=intercept, Boundary.knots =Boundary.knots))/(1e-6)
    } else {
      res <- (splines2::bSpline(t+1e-6, knots=knots, intercept=intercept, Boundary.knots =Boundary.knots)-splines2::bSpline(t-1e-6, knots=knots, intercept=intercept, Boundary.knots =Boundary.knots))/(2*1e-6)
    }
    res
  })
  t(div)
}

# fitLinkLinear impliments the Step 1 of the porposed algorithm.
# stimatedNuissance should contain two lists. One for propensity, the other for the estimated S.
# The predictions nad inputs is always for all the data.
# SplitIndex containes three list. $nuisance, $fit ,$efficient
# treatment is 0 or 1
# covariate does not contain intercept
fitLinkLinear <- function(covariate, response, treatment, estimatedNuisance, splitIndex = NULL, lossType = 'exponential', weights = NULL, tol = 1e-3, intercept = TRUE){
  if(is.null(weights)){
    weights <- rep(1, times=NROW(covariate))
  }
  pi <- estimatedNuisance$pi[splitIndex$fit]
  pi_T <- pi * treatment[splitIndex$fit] + (1-pi) * (1-treatment[splitIndex$fit])
  S <- estimatedNuisance$S[splitIndex$fit]
  trt_train <- 2*(treatment[splitIndex$fit]-0.5)
  x_train <- covariate[splitIndex$fit,]
  y_train <- response[splitIndex$fit]
  w_train <- weights[splitIndex$fit]
  if(intercept){
    x_train <- cbind(1, covariate[splitIndex$fit,])
  }
  obj <- function(beta){
    value <- weighted.mean(y_train/pi_T*loss(trt_train * x_train %*% beta, type = lossType)+trt_train/(2*pi_T)*(S-y_train)*x_train %*% beta, w_train)
    grad <- (apply(apply(x_train,2,function(z){(1/pi_T*y_train*trt_train* loss(trt_train * x_train %*% beta, type = lossType, order = 'derivative')+trt_train/(2*pi_T)*(S-y_train))*z}),
                  2, function(t){weighted.mean(t, w_train)}))
    list(objective = value, gradient = grad)
  }
  beta0 <- array(0,c(NCOL(x_train),1))
  fit <- nloptr::nloptr(x0=beta0, eval_f = obj, opts = list("algorithm"="NLOPT_LD_SLSQP", "xtol_rel"=tol))
  if(intercept){
    fit$solution <- fit$solution[-1]
  }
  fit$solution <- fit$solution/sqrt(sum(fit$solution^2))
  fit <- list(beta=fit$solution, xi=NULL, knods = NULL, solution=NULL)
  fit
}

# updateXi implements the Step 2 of the proposed algorithm
updateXi <- function(fit, covariate, response, treatment, estimatedNuisance, lossType='tanh', splitIndex = NULL, weights = NULL, tol = 1e-3, numberKnots = 5, m_0 =NULL, constraint = TRUE, boundaryPoint=c(-1,1)){
  beta_last <- fit$beta
  if(is.null(weights)){
    weights <- rep(1, times=NROW(covariate))
  }
  pi <- estimatedNuisance$pi[splitIndex$fit]
  pi_T <- pi * treatment[splitIndex$fit] + (1-pi) * (1-treatment[splitIndex$fit])
  S <- estimatedNuisance$S[splitIndex$fit]
  trt_train <- 2*(treatment[splitIndex$fit]-0.5)
  x_train <- covariate[splitIndex$fit,]
  y_train <- response[splitIndex$fit]
  w_train <- weights[splitIndex$fit]
  z_train <- x_train %*% beta_last
  knots <- quantile(z_train, probs=seq(0,1,length.out = numberKnots))
  knots[1] <- knots[1]-0.3*abs(knots[2]-knots[1])
  knots[numberKnots] <- knots[numberKnots]+0.3*abs(knots[numberKnots]-knots[numberKnots-1])
  ns_train <- splines2::bSpline(z_train, knots = knots, intercept = FALSE, Boundary.knots = boundaryPoint)
  ns_train <- cbind(1, ns_train)
  # tuning start
  obj <- function(xi){
    value <- weighted.mean(y_train/pi_T*loss(trt_train * ns_train %*% xi, type = lossType)+trt_train/(2*pi_T)*(S-y_train)*ns_train %*% xi, w_train)
    grad <- (apply(apply(ns_train,2,function(z){(1/pi_T*y_train*trt_train* loss(trt_train * ns_train %*% xi, type = lossType, order = 'derivative')+trt_train/(2*pi_T)*(S-y_train))*z}),
                   2, function(t){weighted.mean(t, w_train)}))
    list(objective = value, gradient = grad)
  }
  #form contrast matrix
  B <- matrix(0,NCOL(ns_train)-1, NCOL(ns_train))
  for(i in 2:(NCOL(ns_train)-1)){
    B[i,i] <- 1
    B[i, i+1] <- -1
  }
  g_ineq <- function(xi){
    value <- c(sum(abs(xi[-1]))-m_0, c(B%*% xi))
    grad <- rbind(c(0,sign(xi[-1])), B)
    list(constraints = value, jacobian = grad)
  }
  if(!constraint){
    g_ineq <- function(xi){
      value <- c(sum(abs(xi[-1]))-m_0)
      grad <- c(0,sign(xi[-1]))
      list(constraints = value, jacobian = grad)
    }
  }
  fit_xi <- nloptr::nloptr(x0=rep(0, times=NCOL(ns_train)), eval_f = obj, eval_g_ineq = g_ineq, opts = list("algorithm"="NLOPT_LD_SLSQP", "xtol_rel"=tol))
  z_test <- seq(min(knots),max(knots), length.out=100)
  ns_test <- splines2::bSpline(z_test, knots = knots, intercept = FALSE, Boundary.knots = boundaryPoint)
  ns_test <- cbind(1, ns_test)
  fit_test <- ns_test %*% fit_xi$solution
  plot(z_test, fit_test)
  fit <- list(beta=beta_last, xi=fit_xi$solution, knots = knots,solution=fit$solution)
  fit
}

# updateBeta implements the Step 3 of the proposed algorithm
updateBeta <- function(fit, covariate, response, treatment, estimatedNuisance, lossType='tanh', splitIndex = NULL, weights = NULL, tol = 1e-3, withConstraint = TRUE, boundaryPoint = c(-1,1), step=0.1){
  beta_last <- fit$beta
  xi <- fit$xi
  if(is.null(weights)){
    weights <- rep(1, times=NROW(covariate))
  }
  pi <- estimatedNuisance$pi[splitIndex$fit]
  pi_T <- pi * treatment[splitIndex$fit] + (1-pi) * (1-treatment[splitIndex$fit])
  S <- estimatedNuisance$S[splitIndex$fit]
  trt_train <- 2*(treatment[splitIndex$fit]-0.5)
  x_train <- covariate[splitIndex$fit,]
  y_train <- response[splitIndex$fit]
  w_train <- weights[splitIndex$fit]
  knots <- fit$knots
  # objective function
  obj <- function(beta){
    z_train <- x_train %*% beta
    ns_train <- splines2::bSpline(z_train, knots = knots, intercept = FALSE, Boundary.knots = boundaryPoint)
    ns_train <- cbind(1, ns_train)
    nsd_train <- cbind(0,nsd(z_train, knots = knots, intercept = FALSE, Boundary.knots = boundaryPoint))
    value <- weighted.mean(y_train/pi_T*loss(trt_train * ns_train %*% xi, type = lossType)+trt_train/(2*pi_T)*(S-y_train)*ns_train %*% xi, w_train)
    grad <- (apply(apply(x_train,2,function(z){(1/pi_T*y_train*trt_train* loss(trt_train * ns_train %*% xi, type = lossType, order = 'derivative')+trt_train/(2*pi_T)*(S-y_train))*nsd_train%*% xi*z}),
                   2, function(t){weighted.mean(t, w_train, na.rm = TRUE)}))
    list(objective = value, gradient = grad)
  }
  constraint <- function(beta){
    return(list('constraints'=sum(beta^2)-1,
           'jacobian'=2*beta))
  }
  if (withConstraint){
    fit_beta <- nloptr::nloptr(x0=fit$beta, eval_f = obj, eval_g_eq = constraint, lb=beta_last-step, ub=beta_last+step, opts = list("algorithm"="NLOPT_LD_SLSQP", "xtol_rel"=tol))
  } else {
    fit_beta <- nloptr::nloptr(x0=fit$beta, eval_f = obj, opts = list("algorithm"="NLOPT_LD_SLSQP", "xtol_rel"=tol))
  }
  fit <- list(beta=fit_beta$solution, xi=fit$xi, knots = fit$knots)
  fit
}

