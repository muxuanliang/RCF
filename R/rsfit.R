# rsfitSolver obtains the optimizer of the propoed optimization given M_n
rsfitSolver <- function(covariate, response, treatment, estimatedNuisance, splitIndex = NULL, lossType = 'logistic', weights = NULL, tol = 10^(-3), m_0=0, constraint = TRUE, boundaryPoint=c(-1,1), numberKnots = 5){
  # initialization
  fit_earl <- ITRFitAll(list(predictor = covariate, treatment = (treatment > 0), outcome = response), loss = 'logistic', sampleSplitIndex = splitIndex$nuisance, propensity = estimatedNuisance$pi, outcome = estimatedNuisance, test=FALSE)
  #fit_init <- fitLinkLinear(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, lossType = lossType, weights = weights, tol = tol)
  iter <- 0

  # iteration start
  fit_last <- NULL
  tmp_beta<- (fit_earl$fit[[1]]$fit$beta[,fit_earl$fit[[1]]$fit$lambda==fit_earl$fit[[1]]$fit$lambda.min] + fit_earl$fit[[2]]$fit$beta[,fit_earl$fit[[2]]$fit$lambda==fit_earl$fit[[2]]$fit$lambda.min])/2
  if (sum(tmp_beta^2)<10^-3){
    fit_last$beta <- (tmp_beta+0.1)/sqrt(sum((tmp_beta+0.1)^2))
  } else {
    fit_last$beta <- tmp_beta/sqrt(sum(tmp_beta^2))
  }

  fit_last$xi <- 0
  diff <- 1
  step <- 0.1
  while((diff > tol) & (iter < 100)){
    if(any(is.na(fit_last$beta))){
      print(fit_last)
      print(fit_xi)
      print(fit_beta)
      print(iter)
    }
    # updateXi
    fit_xi <- updateXi(fit_last, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = 1e-5, m_0=m_0, constraint=constraint, boundaryPoint = boundaryPoint, numberKnots = numberKnots)
    # updateBeta
    fit_beta <- updateBeta(fit_xi, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol, withConstraint = TRUE, boundaryPoint = boundaryPoint, step=step)
    # update no constraint
    #fit_tmp <- updateBeta(fit_xi, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol, withConstraint = FALSE, boundaryPoint = boundaryPoint)
    #update

    diff_coef <- max(abs(fit_beta$beta-fit_last$beta))
    diff_xi <- max(c(abs(fit_beta$xi-fit_last$xi))) #abs(fit_beta$beta-fit_tmp$beta))
    diff <- max(diff_coef, diff_xi)
    iter <- iter+1
    fit_last <- fit_beta
  }
  # final update
  fit <- updateXi(fit_last, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol, m_0=m_0, constraint=constraint, boundaryPoint = boundaryPoint, numberKnots = numberKnots)
  fit_tmp <- updateBeta(fit, covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex,lossType = lossType, weights = weights, tol = tol, withConstraint = FALSE, boundaryPoint = boundaryPoint)
  fit$solution <- fit_tmp$beta
  fit$boundaryPoint <- boundaryPoint
  fit
}

# rsfitSplit obtain the estimation of RSICF assuming a single index model on a slited data. The variance estimaor of the coefficiencts are provided.
rsfitSplit <- function(covariate, response, treatment, splitIndex = NULL, propensityModel = 'glmnet', estimatedPropensity = NULL, outcomeModel = 'kernel', estimatedOutcome = NULL, lossType = 'logistic', weights = NULL, tol = 1e-3, propensityFormula=NULL, outcomeFormula=NULL, constraint=TRUE, boundaryPoint=c(-1,1)){
  # fit nuisance parameter
  estimatedNuisance <- NULL
  estimatedNuisance$pi <- estimatedPropensity
  if (is.null(estimatedPropensity)){
    data <- list(predictor = covariate, treatment = treatment)
    predictedPropensityAll <- getPropensityModel(data, method = propensityModel, splitIndex = splitIndex, Formula = propensityFormula, predictAll = TRUE)
    estimatedNuisance$pi <- predictedPropensityAll
  }
  estimatedNuisance$S <- estimatedOutcome$control+estimatedOutcome$treatment
  estimatedNuisance$control <- estimatedOutcome$control
  estimatedNuisance$treatment <- estimatedOutcome$treatment
  if (is.null(estimatedOutcome)){
    data <- list(predictor = covariate, treatment = treatment, outcome = response)
    predictedOutcomeAll <- getOutcomeModel(data, method = outcomeModel, splitIndex = splitIndex, Formula = outcomeFormula, predictAll = TRUE)
    estimatedNuisance$S <- predictedOutcomeAll$control+predictedOutcomeAll$treatment
    estimatedNuisance$control <- predictedOutcomeAll$control
    estimatedNuisance$treatment <- predictedOutcomeAll$treatment
  }

  # fit conditional variance
  residual_square <- (response - (predictedOutcomeAll$control * (treatment==0) + predictedOutcomeAll$treatment * (treatment==1)))^2
  data <- list(predictor = covariate, treatment = treatment, outcome = residual_square)
  predictedVarAll <- getOutcomeModel(data, method = outcomeModel, splitIndex = splitIndex, Formula = outcomeFormula, predictAll = TRUE)
  estimatedNuisance$var_control <- predictedVarAll$control
  estimatedNuisance$var_treated <- predictedVarAll$treatment

  # weights
  if(is.null(weights)){
    weights <- rep(1, times=NROW(covariate))
  }

  # solve
  m_0 <- 5
  numberKnots <- 3
  fit_hat <- rsfitSolver(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, lossType = lossType, weights =weights, tol = tol, m_0=m_0, constraint = constraint, boundaryPoint = boundaryPoint, numberKnots = numberKnots)
  diff_min <- max(abs(fit_hat$solution-fit_hat$beta))
  fit_hat_min <- fit_hat
  dif_tmp <- diff_min
  while((diff_min>10^(-3)) & (numberKnots <= 7) & (m_0 <= 2^5)){
    early_step <- FALSE
    if (sum(abs(fit_hat$xi[-1]))<m_0){
      numberKnots <- numberKnots + 1
      early_step <- TRUE
    } else {
      m_0 <- 2*m_0
    }
    fit_hat <- rsfitSolver(covariate, response, treatment, estimatedNuisance, splitIndex = splitIndex, lossType = lossType, weights =weights, tol = tol, m_0=m_0, constraint = constraint, boundaryPoint = boundaryPoint, numberKnots = numberKnots)
    dif_tmp <- max(abs(fit_hat$solution-fit_hat$beta))
    if (dif_tmp < diff_min){
      fit_hat_min <- fit_hat
      diff_min <- dif_tmp
    }
    if ((dif_tmp > diff_min)&early_step&(sum(abs(fit_hat$xi[-1]))<m_0)){
      break
    }
    if ((dif_tmp > diff_min)&(!early_step)&(sum(abs(fit_hat$xi[-1]))>=m_0)){
      break
    }
  }

  # return
  list(fit=fit_hat_min, covariate=covariate, response=response, treatment=treatment, estimatedNuisance=estimatedNuisance, splitIndex=splitIndex, lossType = lossType, weights =weights)
}

# rsfit obtains sysmetric results
rsfit <- function(covariate, response, treatment, splitIndex = NULL, propensityModel = 'glmnet', estimatedPropensity = NULL, outcomeModel = 'kernel', estimatedOutcome = NULL, lossType = 'logistic', weights = NULL, tol = 1e-3, propensityFormula=NULL, outcomeFormula=NULL, parallel = FALSE, constraint = TRUE, boundaryPoint = c(-1,1), efficient = TRUE, local=TRUE){
  n <- NROW(covariate)
  if(is.null(splitIndex)){
    random_index <- sample(c(1,2,3), n, replace = TRUE)
    splitIndex <- NULL
    splitIndex$nuisance <- (random_index == 1)
    splitIndex$fit <- (random_index == 2)
    splitIndex$infer <- (random_index == 3)
  }
  fit <- NULL
  fit$fit <- NULL
  if(!parallel){
    fit_tmp <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
    fit$fit[[1]] <- rsInference(fit_tmp, efficient=efficient, local = local)
    # rotate the index
    tmp <- splitIndex$nuisance
    splitIndex$nuisance <- splitIndex$infer
    splitIndex$infer <- splitIndex$fit
    splitIndex$fit <- tmp

    fit_tmp <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
    fit$fit[[2]] <- rsInference(fit_tmp, efficient=efficient, local = local)
    # rotate the index
    tmp <- splitIndex$nuisance
    splitIndex$nuisance <- splitIndex$infer
    splitIndex$infer <- splitIndex$fit
    splitIndex$fit <- tmp

    fit_tmp <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
    fit$fit[[3]] <- rsInference(fit_tmp, efficient=efficient, local = local)
  } else {
    library(doParallel)
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
    cl <- makeCluster(3)
    registerDoParallel(cl)
    result <- foreach(iter = c(1,2,3), .packages = c("RSICF","ITRInference", 'glmnet', 'mgcv'))%dopar%{
      splitIndex_local <- NULL
      if (iter == 1){
        splitIndex_local <- splitIndex
      } else if (iter == 2) {
        splitIndex_local$nuisance <- splitIndex$infer
        splitIndex_local$infer <- splitIndex$fit
        splitIndex_local$fit <- splitIndex$nuisance
      } else {
        splitIndex_local$nuisance <- splitIndex$fit
        splitIndex_local$infer <- splitIndex$nuisance
        splitIndex_local$fit <- splitIndex$infer
      }
      fit_tmp <- rsfitSplit(covariate, response, treatment, splitIndex = splitIndex_local, lossType = lossType, propensityModel = propensityModel, estimatedPropensity = estimatedPropensity, outcomeModel=outcomeModel, estimatedOutcome = estimatedOutcome, weights = weights, tol=tol, propensityFormula = propensityFormula, outcomeFormula = outcomeFormula, constraint = constraint, boundaryPoint = boundaryPoint)
      rsInference(fit_tmp, efficient=efficient, local = local)
    }
    stopCluster(cl)
    fit$fit <- result
  }
  fit$betaAN <- (fit$fit[[1]]$betaAN+fit$fit[[2]]$betaAN+fit$fit[[3]]$betaAN)/3
  fit$sigmaAN <- (fit$fit[[1]]$sigmaAN+fit$fit[[2]]$sigmaAN+fit$fit[[3]]$sigmaAN)/3

  # normalize
  fit$betaAN.renorm <- (fit$fit[[1]]$betaAN.renorm+fit$fit[[2]]$betaAN.renorm+fit$fit[[3]]$betaAN.renorm)/3
  fit$sigmaAN.renorm <- (fit$fit[[1]]$sigmaAN.renorm+fit$fit[[2]]$sigmaAN.renorm+fit$fit[[3]]$sigmaAN.renorm)/3
  norm1 <- sqrt(sum(fit$betaAN.renorm^2))
  trans <- matrix(0, NROW(fit$sigmaAN.renorm), NCOL(fit$sigmaAN.renorm))
  for (i in 1: NROW(fit$sigmaAN.renorm)){
    for (j in 1:NCOL(fit$sigmaAN.renorm)){
      if (i==j){
        trans[i,j] <- (norm1^2-(fit$betaAN.renorm[i])^2)/(norm1^3)
      } else {
        trans[i,j] <- -(fit$betaAN.renorm[i]) * (fit$betaAN.renorm[j])/(norm1^3)
      }
    }
  }
  fit$sigmaAN.renorm <- trans %*% fit$sigmaAN.renorm %*% trans
  fit$betaAN.renorm <- fit$betaAN.renorm/norm1

  #return
  fit$sigmaAN <- diag(fit$sigmaAN)
  fit$sigmaAN.renorm <- diag(fit$sigmaAN.renorm)
  fit
}





