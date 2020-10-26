rsInference <- function(rsfit, efficient = FALSE, local=TRUE){
  # set parameter
  estimatedNuisance <- rsfit$estimatedNuisance
  splitIndex <- rsfit$splitIndex
  covariate <- rsfit$covariate
  response <- rsfit$response
  treatment <- rsfit$treatment
  lossType <- rsfit$lossType

  # set nuisance parameter
  pi <- estimatedNuisance$pi
  S <- estimatedNuisance$S

  # set weight
  weights <- 1
  if(efficient & local){
    eta <- predict.rsfitSplit(rsfit$fit, newx = covariate, type='eta')
    theta <- predict.rsfitSplit(rsfit$fit, newx = covariate)
    # augmentation
    ratio <- (1-eta)/(1+eta)
    inv_W <- (1/estimatedNuisance$pi * ratio + 1/(1-estimatedNuisance$pi) * 1/ratio)*(exp(theta/2)+exp(-theta/2))^2
    weights <- inv_W^(-1) * estimatedNuisance$S
  }
  if(efficient & (!local)){
    eta <- predict.rsfitSplit(rsfit$fit, newx = covariate, type='eta')
    theta <- predict.rsfitSplit(rsfit$fit, newx = covariate)
    # augmentation
    ratio <- (1-eta)/(1+eta)
    var_treated <- estimatedNuisance$var_treated
    var_control <- estimatedNuisance$var_control
    inv_W <- (1/estimatedNuisance$pi * ratio * var_treated + 1/(1-estimatedNuisance$pi) * 1/ratio * var_control)*(exp(theta/2)+exp(-theta/2))^2
    weights <- inv_W^(-1) * estimatedNuisance$S
  }

  # set bias of x
  link <- covariate %*% rsfit$fit$beta
  exp_gam <- function(x, y){
    shadow_data <- data.frame(predictor = x, response = y)
    shadow_model <- gam(response~predictor, data = shadow_data)
    as.matrix(predict(shadow_model, newdata = data.frame(predictor = link),type='response'))
  }
  tmp_w <- response/(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment))
  ESX <- apply(covariate[splitIndex$fit,], 2, function(t){exp_gam(link[splitIndex$fit], weights[splitIndex$fit]*tmp_w[splitIndex$fit] * t)})
  ES <- exp_gam(link[splitIndex$fit], tmp_w[splitIndex$fit]*weights[splitIndex$fit])

  # start inference
  predict <- predict.rsfitSplit(rsfit$fit, newx = covariate)
  tildeX <- apply(ESX,2,function(t){t/ES})
  w0 <- predict.rsfitSplit(rsfit$fit, newx = covariate, derivative=TRUE)
  w1 <- response/(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment)) * loss(2*(treatment-0.5) * predict, type = lossType, order = 'hessian') * w0^2 * weights
  w2 <- (2*(treatment-0.5)*response/(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment)) * loss(2*(treatment-0.5) * predict, type = lossType, order = 'derivative')+2*(treatment-0.5)*(estimatedNuisance$S-response)/(2*(estimatedNuisance$pi*treatment+(1-estimatedNuisance$pi)*(1-treatment)))) * w0 * weights
  sampleIndex <- c(1:NROW(covariate))
  W1 <- matrix(0, NCOL(covariate), NCOL(covariate))
  W2 <- matrix(0, NCOL(covariate), NCOL(covariate))
  S <- matrix(0, NCOL(covariate), 1)
  for (iter in sampleIndex[splitIndex$fit]){
    W1 <- W1+w1[iter] * (covariate[iter,]) %*% t(covariate[iter,]-tildeX[iter,])
    W2 <- W2+(w2[iter])^2 * (covariate[iter,]-tildeX[iter,]) %*% t(covariate[iter,]-tildeX[iter,])
    S <- S+w2[iter]* (covariate[iter,]-tildeX[iter,])
  }
  W1 <- W1/sum(splitIndex$fit)
  #W1 <- W1 + 10^(-5)*diag(min(abs(diag(W1))),NCOL(covariate),NCOL(covariate))
  W2 <- W2/sum(splitIndex$fit)
  S <- S/sum(splitIndex$fit)

  # deal with ill-posed variance
  t <- try(solve(W1), silent = TRUE)
  while('try-error' %in% class(t)){
    W1 <- W1 + 1e-5 * diag(1,NCOL(covariate),NCOL(covariate))
    t <- try(solve(W1), silent = TRUE)
  }

  betaAN <- rsfit$fit$beta-solve(W1) %*% S
  sigmaAN <- solve(W1) %*% W2 %*% solve(W1)

  # if sigmaAN < 0
  if (min(diag(sigmaAN))<0){
    W1 <- matrix(0, NCOL(covariate), NCOL(covariate))
    for (iter in sampleIndex[splitIndex$fit]){
      W1 <- W1+w1[iter] * (covariate[iter,]-tildeX[iter,]) %*% t(covariate[iter,]-tildeX[iter,])
    }
    W1 <- W1/sum(splitIndex$fit)
    t <- try(solve(W1), silent = TRUE)
    while('try-error' %in% class(t)){
      W1 <- W1 + 1e-5 * diag(1,NCOL(covariate),NCOL(covariate))
      t <- try(solve(W1), silent = TRUE)
    }
    betaAN <- rsfit$fit$beta-solve(W1) %*% S
    sigmaAN <- solve(W1) %*% W2 %*% solve(W1)
    print('small variance')
  }

  # normalize
  norm1 <- sqrt(sum(betaAN^2))
  trans <- matrix(0, NROW(sigmaAN), NCOL(sigmaAN))
  for (i in 1: NROW(sigmaAN)){
    for (j in 1:NCOL(sigmaAN)){
      if (i==j){
        trans[i,j] <- (norm1^2-(betaAN[i])^2)/norm1^3
      } else {
        trans[i,j] <- -(betaAN[i]) * (betaAN[j])/norm1^3
      }
    }
  }
  sigmaAN.renorm <- trans %*% sigmaAN %*% trans
  betaAN.renorm <- betaAN/norm1
  #return
  list(fit=rsfit$fit, betaAN=betaAN, W1=W1, W2=W2, sigmaAN=sigmaAN, betaAN.renorm=betaAN.renorm, sigmaAN.renorm=sigmaAN.renorm)
}
