
prepare_output <- function(start_time, estimates, obs, mask, beta=NA, beta.estim=NA, M=NA, M.estim=NA, 
                           LogLik_SI=NA, test_error = utils$error_metric$rmse){
  estim.test <- estimates[mask==0]
  estim.train <- estimates[mask != 0]
  obs.train <- obs[mask!=0]
  obs.test <- obs[mask==0]
  list(
    time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs"))),
    error.test = test_error(estim.test, obs.test),
    corr.test = cor(estim.test, obs.test),
    error.train = test_error(estim.train, obs.train),
    error.M = tryCatch(test_error(M.estim, M), error = function(x) NA),
    error.beta = tryCatch(test_error(beta.estim, beta), error = function(x) NA),
    rank_M = tryCatch(qr(M.estim)$rank, error = function(x) NA),
    rank_beta = tryCatch(qr(beta.estim)$rank, error = function(x) NA),
    
    sparse_in_sparse = tryCatch(sum(beta == 0 & beta.estim == 0) /
                                  (sum(beta == 0) +  1e-17), error = function(x) NA),
    nonsparse_in_nonsparse = tryCatch(sum(beta != 0 & beta.estim != 0) /
                                        (sum(beta != 0) +  1e-17), error = function(x) NA)
  ) -> results
  
  
  if(is.null(LogLik_SI)){
    results$likelihood_ratio_index <- NA
    results$Cox_Snell_R2 <- NA
  }else{
    residuals <- obs.test - estim.test
    LogLik <- utils$logLikelihood(residuals)
    n <- length(residuals)
    results$likelihood_ratio_index <- utils$Likelihood_ratio_index(LogLik, LogLik_SI)
    results$Cox_Snell_R2 <-utils$ Cox_Snell_R2(LogLik, LogLik_SI, n)
  }
  return(results)
}




Mao_Sim_Wrapper <-
  function(dat,
           lambda.1_grid = seq(0, 1, length = 20),
           lambda.2_grid = seq(0.9, 0.1, length = 20),
           alpha_grid = c(1),
           ncores = 1,
           # keep it > 1
           n_folds = 5,
           weight_function = Mao_weights$uniform,
           LogLik_SI = NULL,
           ...) {
    start_time = Sys.time()
    fiti <- Mao.cv(
      Y = dat$Y,
      X = dat$X,
      W = dat$W,
      n_folds = n_folds,
      lambda.1_grid = lambda.1_grid,
      lambda.2_grid = lambda.2_grid,
      alpha_grid = alpha_grid,
      seed = 2023,
      numCores = ncores,
      n1n2_optimized = TRUE,
      test_error = utils$error_metric$rmse,
      theta_estimator = weight_function,
      sequential = FALSE
    )
    
    fit. <- fiti$fit
    results = list(model = "Mao")
    results$lambda.beta = fiti$best_parameters$lambda.1
    results$lambda.M = fiti$best_parameters$lambda.2
    results <- c(results,
                 prepare_output(start_time, fit.$estimates, dat$O, dat$W, dat$beta, fit.$beta, dat$M, fit.$M, LogLik_SI))    
    
    
    results
}




SImpute_Sim_Wrapper <- function(dat, ...) {
  start_time = Sys.time()
  fit. <- simpute.cv(
    Y_train = as.matrix(dat$fit_data$train),
    y_valid = dat$fit_data$valid,
    W_valid = dat$fit_data$W_valid,
    y = dat$Y,
    n.lambda = 20,
    trace = FALSE,
    print.best = FALSE,
    tol = 5,
    thresh = 1e-6,
    rank.init = 2,
    rank.limit = 30,
    rank.step = 2,
    maxit = 600,
    seed = NULL
  )
  results = list(model = "SoftImpute")
  results$lambda.beta = NA
  results$lambda.M = fit.$lambda
  results <- c(results,
                        prepare_output(start_time, fit.$estimates, dat$O, dat$W, M.estim = fit.$estimates))    
  
  LogLik <- utils$logLikelihood(dat$O[dat$W==0] - fit.$estimates[dat$W==0])
  return(list(results=results, LogLik=LogLik))
}

#--------------------------------------------------------------------------------------
CASMC_Ridge_Sim_Wrapper <-
  function(dat,
           max_cores = 20,
           LogLik_SI = NULL,
           ...) {
    start_time = Sys.time()
    
    fiti <- CASMC_Ridge_cv(
      y_train = dat$fit_data$train,
      X = dat$X,
      y_valid = dat$fit_data$valid,
      W_valid = dat$fit_data$W_valid,
      y = dat$fit_data$Y,
      error_function = utils$error_metric$rmse,
      thresh = 1e-6,
      maxit = 300,
      trace = FALSE,
      print.best = F,
      quiet = FALSE,
      warm = NULL,
      track = F,
      max_cores = max_cores,
      seed = NULL
    )
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC-Ridge")
    results$lambda.beta = fiti$lambda.beta
    results$lambda.M = fit.$lambda
    results <- c(results,
                          prepare_output(start_time, fit.$estimates, dat$O, dat$W, dat$beta, fit.$beta, dat$M, fit.$M, LogLik_SI))  
    results
  }
#-------


CASMC_Nuclear_Sim_Wrapper <-
  function(dat,
           LogLik_SI = NULL,
           ...) {
    start_time = Sys.time()
    
    fiti <- CASMC_Nuclear_cv(
      y_train = dat$fit_data$train,
      X = dat$X,
      y_valid = dat$fit_data$valid,
      W_valid = dat$fit_data$W_valid,
      y = dat$fit_data$Y,
      error_function = utils$error_metric$rmse,
      warm = NULL,
      trace = F,
      quiet = T,
      track = F,
      step3 = T,
      use_warmstart = T,
      seed = NULL
    )
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = utils$unsvd(fit.)
    fit.$beta = utils$unsvd(fit.$beta)
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC-Nuclear")
    results$lambda.beta = fiti$hparams$lambda.beta
    results$lambda.M = fiti$hparams$lambda.M
    
    results <- c(results,
                          prepare_output(start_time, fit.$estimates, dat$O, dat$W, dat$beta, fit.$beta, dat$M, fit.$M, LogLik_SI))  
    
    results
  }
#----------------------------------------------------
CASMC_Lasso_Sim_Wrapper <-
  function(dat,
           max_cores = 20,
           LogLik_SI = NULL,
           ...) {
    start_time = Sys.time()
    fiti <- CASMC_Lasso_cv(
      y_train = dat$fit_data$train,
      X = dat$X,
      y_valid = dat$fit_data$valid,
      W_valid = dat$fit_data$W_valid,
      y = dat$fit_data$Y,
      trace = 0,
      print.best = F,
      warm = NULL,
      quiet = T, 
      max_cores = max_cores
    ) 
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC-Lasso")
    results$lambda.beta = fiti$hparams$lambda.beta
    results$lambda.M = fiti$hparams$lambda.M
    
    results <- c(results,
                          prepare_output(start_time, fit.$estimates, dat$O, dat$W, dat$beta, fit.$beta, dat$M, fit.$M, LogLik_SI))  
    
    results
  }
#----------------------------------------------
Naive_Sim_Wrapper <- function(dat, ...) {
  start_time = Sys.time()
  fit. <- naive_fit(dat$Y, dat$X)
  results = list(model = "Naive")
  results$lambda.beta = NA
  results$lambda.M = NA
  results <- c(results,
                        prepare_output(start_time, fit.$estimates, dat$O, dat$W, dat$beta, fit.$beta, dat$M, fit.$M))  
  results
}
