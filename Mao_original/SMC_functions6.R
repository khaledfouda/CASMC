
cv_validate_sim = function(diagD, tau2_grid, tau1_grid_length, svdY, alpha, uvalidatemat, 
                           Xbeta_validate, umat_validate) {
  
  # This function calculates the validation error to tuning parameters.
  
  cv_validate = array(0, c(length(tau2_grid), tau1_grid_length))
  
  for (i in 1:length(tau2_grid)) {
    
    B_validate = as.vector(SVTE_alpha(svdY$u, svdY$d, svdY$v, tau2_grid[i], 
                                      alpha))[uvalidatemat != 0]
    
    for (j in 1:tau1_grid_length) {
      
      cv_validate[i, j] = sum((Xbeta_validate[j, ] + B_validate - umat_validate)^2)/sqrt(sum(uvalidatemat != 
                                                                                               0))
      
    }
  }
  
  return(cv_validate)
}


########################## Scale SVT
SVTE_alpha = function(svdu, svdd, svdv, tau_svd_ratio, alpha_ratio) {
  
  # This function performs the singular value soft-thresholding and scaling
  # procedures.
  
  return(svdu %*% (pmax(svdd - svdd[1] * tau_svd_ratio * alpha_ratio, 0) * t(svdv))/(1 + 
                                                                                       2 * svdd[1] * tau_svd_ratio * (1 - alpha_ratio)))
  
}


SMCfit = function(Aobs, X, tau_beta_ratio, tau_svd_ratio, alpha_ratio) {
  
  #' This function performs the matrix completion part in proposed method with fix tunning parameters.
  #'
  #' Arguments:
  #'    Aobs:             The observation sample matrix.
  #'    X:                The sample covariate matrix.
  #'    diagD:            The input diagonal matrix with inclusion probabilities pi_{i}.
  #'    tau_beta_ratio:   The input tunning parameter 'tau_{1}'.
  #'    tau_svd_ratio:    The input tunning parameter 'tau_{2}'.
  #'    alpha_ratio:      The tunning parameter 'alpha'.
  #'
  #' Outputs: 
  #'    A list is returned, with values: Estimated Ahat, Bhat, betahat and rank.
  
  n1 = dim(Aobs)[1]
  n2 = dim(Aobs)[2]
  m = dim(X)[2]
  
  iomega = as.numeric(Aobs != 0)
  omega = matrix(iomega, n1, n2) # W
  #---------------------------------
  # estimating Theta (named gammahat)
  #------------------------------------
  Xtheta = X#[, 1:3]# ???
  
  Xadj = X
  PX = Xadj %*% solve(t(Xadj) %*% Xadj) %*% t(Xadj)
  Eye = diag(1, n1)
  PXp = Eye - PX
  
  # the whole gammahat ???
  gammahat = matrix(rep(NA, n1 * n2), nrow = n1)
  for (i in 1:n2) { # is this by columns? should we only conser column i in X?
    
    thetadata = data.frame(cbind(omega[, i], Xtheta)) # ???
    thetaglmout = glm(X1 ~ ., family = binomial(logit), data = thetadata)
    gammahat[, i] = predict(thetaglmout, type = "response")
  }
  #-----------------------------------------------------------
  
  Aobsw = (Aobs/gammahat)
  
  PXpAobsw = PXp %*% Aobsw
  
  ### The following part is different from the original one.
  betahat = solve(t(Xadj) %*% Xadj + svd(t(Xadj) %*% Xadj)$d[1] * tau_beta_ratio * 
                    diag(1, m)) %*% t(Xadj) %*% Aobsw
  Xbetahat = X %*% betahat
  
  svdPXpAobsw = svd(PXpAobsw)
  
  Bhat = SVTE_alpha(svdPXpAobsw$u, svdPXpAobsw$d, svdPXpAobsw$v, tau_svd_ratio, 
                    alpha_ratio)
  
  rankB = sum(pmax(svdPXpAobsw$d - svdPXpAobsw$d[1] * tau_svd_ratio * alpha_ratio, 
                   0) > 0)
  
  ###################### 
  Ahat = Xbetahat + Bhat
  
  return(list(Ahat = Ahat, Bhat = Bhat, betahat = betahat, rank = rankB + m))
}


SMCfit_cv = function(Y, X.cov, mask, Aobs, nfolds = 5, 
                     tau1_grid = seq(0, 1, length = 30), 
                     tau2_grid = seq(0.9, 0.1, length = 30), alpha_grid = seq(0.992, 1, length = 20),seed=2023) {
  
  #' This function performs the matrix completion part in the proposed method with the tuning parameter chosen by cross-validation.
  #'
  #' Arguments:
  #'    Y:                The population matrix containing the values of interest.
  #'    X.cov:            The population matrix of covariate.
  #'    prob.mis:         The missing probability for each element of the population.
  #'    pps.prob:         The normalized inclusion probability for probability-proportional-to-size sampling.
  #'    nfolds:           The number of cross-validation folds; with default value 5.
  #'    tau1_grid:        The grid used in search for parameter '$\tau_{1}$'; with default value seq(0,1,length=30).
  #'    tau2_grid:        The grid used in search for parameter '$\tau_{2}$'; with default value seq(0.9,0.1,length=30).
  #'    alpha_grid:       The grid used for alpha; with default value seq(0.992,1,length=20).
  #'
  #' Outputs: 
  #'    A list is returned, with values: Estimated Ahat, Bhat, betahat and rank.
  
  set.seed(seed)
  N.popu = nrow(Y)
  
  
  
  sample.full = Y
  x.sample = X.cov
  
  omega = mask
  
  
  #Aobs = Y * omega
  n1 = dim(Aobs)[1]
  n2 = dim(Aobs)[2]
  m = dim(x.sample)[2]
  
  Xtheta = x.sample#[, 1:3]
  
  Xadj = x.sample
  PX = Xadj %*% solve(t(Xadj) %*% Xadj) %*% t(Xadj)
  Eye = diag(1, n1)
  PXp = Eye - PX
  
  iomega = as.numeric(omega != 0)
  iobs = which(iomega == 1)
  imiss = which(iomega == 0)
  
  gammahat = matrix(rep(NA, n1 * n2), nrow = n1)
  for (i in 1:n2) {
    
    thetadata = data.frame(cbind(omega[, i], Xtheta))
    thetaglmout = glm(X1 ~ ., family = binomial(logit), data = thetadata)
    gammahat[, i] = predict(thetaglmout, type = "response")
    
  }
  
  Aobsw = (Aobs/gammahat)
  
  PXpAobsw = PXp %*% Aobsw
  
  svdPXpAobsw = svd(PXpAobsw)
  
  folds <- cut(sample(seq(1, length(iobs))), breaks = nfolds, labels = FALSE)
  iomegatrain = matrix(0, n1 * n2, nfolds)
  iomegavalidate = matrix(0, n1 * n2, nfolds)
  omegatrain = array(NA, c(n1, n2, nfolds))
  omegavalidate = array(NA, c(n1, n2, nfolds))
  
  gammatrainhat = array(NA, c(n1, n2, nfolds))
  Aobswtrain = array(NA, c(n1, n2, nfolds))
  PXpAobswtrain = array(NA, c(n1, n2, nfolds))
  svdPXpAobswtrain = list()
  A_validate = list()
  Xbetavalidate = list()
  
  for (fold_num in 1:nfolds) {
    
    validateIndexes <- which(folds == fold_num, arr.ind = TRUE)
    ivalidate <- iobs[validateIndexes]
    itrain <- iobs[-validateIndexes]
    
    iomegatrain[itrain, fold_num] = 1
    omegatrain[, , fold_num] = matrix(iomegatrain[, fold_num], n1, n2)
    
    iomegavalidate[ivalidate, fold_num] = 1
    omegavalidate[, , fold_num] = matrix(iomegavalidate[, fold_num], n1, n2)
    
    for (i in 1:n2) {
      
      thetadata = data.frame(cbind(omegatrain[, i, fold_num], Xtheta))
      thetaglmout = glm(X1 ~ ., family = binomial(logit), data = thetadata)
      gammatrainhat[, i, fold_num] = predict(thetaglmout, type = "response")
      
    }
    
    Aobswtrain[, , fold_num] = omegatrain[, , fold_num] * Aobs/gammatrainhat[, 
                                                                             , fold_num]
    PXpAobswtrain[, , fold_num] = PXp %*% Aobswtrain[, , fold_num]
    svdPXpAobswtrain = append(svdPXpAobswtrain, list(svd(PXpAobswtrain[, , fold_num])))
    
    A_validate = append(A_validate, list(as.vector(Aobs)[iomegavalidate[, fold_num] != 
                                                           0]))
    
    Xbetavalidate = append(Xbetavalidate, list(matrix(nrow = length(tau1_grid), 
                                                      ncol = sum(iomegavalidate[, fold_num]))))
  }
  
  for (j in 1:length(tau1_grid)) {
    
    PX = Xadj %*% solve(t(Xadj) %*% Xadj + svd(t(Xadj) %*% Xadj)$d[1] * tau1_grid[j] * 
                          diag(1, m)) %*% t(Xadj)
    
    for (fold_num in 1:nfolds) {
      Xbetavalidate[[fold_num]][j, ] = as.vector(PX %*% Aobswtrain[, , fold_num])[iomegavalidate[, 
                                                                                                 fold_num] != 0]
    }
  }
  
  cv = array(0, c(nfolds, length(tau2_grid), length(tau1_grid), length(alpha_grid)))
  for (fold_num in 1:nfolds) {
    for (alpha_num in 1:length(alpha_grid)) {
      
      cv[fold_num, , , alpha_num] = cv_validate_sim(diagD, tau2_grid, length(tau1_grid), 
                                                    svdPXpAobswtrain[[fold_num]],
                                                    alpha_grid[alpha_num], iomegavalidate[, fold_num],
                                                    Xbetavalidate[[fold_num]], A_validate[[fold_num]])
      
    }
  }
  
  cv_grid = which(apply(cv, 2:4, sum) == min(apply(cv, 2:4, sum)), arr.ind = T)[1, 
  ]
  
  return(list(cv = cv, Aobs = Aobs, cv_min = min(apply(cv, 2:4, sum)), tau_beta_ratio = tau1_grid[cv_grid[2]], tau_svd_ratio = tau2_grid[cv_grid[1]], 
              alpha_ratio = alpha_grid[cv_grid[3]]))
}
