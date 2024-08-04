
prepare_output <- function(start_time, estimates, obs, mask, beta=NA, beta.estim=NA, M=NA, M.estim=NA, LogLik_SI=NA, test_error = error_metric$rmse){
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
    LogLik <- logLikelihood(residuals)
    n <- length(residuals)
    results$likelihood_ratio_index <- Likelihood_ratio_index(LogLik, LogLik_SI)
    results$Cox_Snell_R2 <- Cox_Snell_R2(LogLik, LogLik_SI, n)
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
      test_error = error_metric$rmse,
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
  
  LogLik <- logLikelihood(dat$O[dat$W==0] - fit.$estimates[dat$W==0])
  return(list(results=results, LogLik=LogLik))
}

#--------------------------------------------------------------------------------------
CASMC_0_Sim_Wrapper <-
  function(dat,
           max_cores = 20,
           LogLik_SI = NULL,
           ...) {
    start_time = Sys.time()
    
    fiti <- CASMC0_cv(
      y_train = dat$fit_data$train,
      X = dat$X,
      y_valid = dat$fit_data$valid,
      W_valid = dat$fit_data$W_valid,
      y = dat$fit_data$Y,
      error_function = error_metric$rmse,
      lambda.factor = 1 / 4,
      lambda.init = NULL,
      n.lambda = 20,
      rank.init = 2,
      rank.limit = 30,
      rank.step = 2,
      pct = 0.98,
      lambda.a = 0,
      S.a = NULL,
      lambda.b = 0,
      S.b = NULL,
      early.stopping = 1,
      thresh = 1e-6,
      maxit = 300,
      trace = FALSE,
      print.best = F,
      quiet = FALSE,
      warm = NULL,
      lambda.beta.grid = "default",
      track = F,
      max_cores = max_cores,
      seed = NULL
    )
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC-0")
    results$lambda.beta = fiti$lambda.beta
    results$lambda.M = fit.$lambda
    results <- c(results,
                          prepare_output(start_time, fit.$estimates, dat$O, dat$W, dat$beta, fit.$beta, dat$M, fit.$M, LogLik_SI))  
    results
  }
#-------


CASMC_2_Sim_Wrapper <-
  function(dat,
           LogLik_SI = NULL,
           ...) {
    start_time = Sys.time()
    
    fiti <- CASMC2_cv2(
      y_train = dat$fit_data$train,
      X = dat$X,
      y_valid = dat$fit_data$valid,
      W_valid = dat$fit_data$W_valid,
      y = dat$fit_data$Y,
      error_function = error_metric$rmse,
      warm = NULL,
      M_cv_param = list(
        rank.init = 2,
        rank.limit = 30,
        rank.step = 2,
        pct = 0.98,
        lambda.factor = 1/4,
        lambda.init = NULL,
        n.lambda = 20, 
        early.stopping = 1
      ),
      beta_cv_param = list(
        rank.init = 2,
        rank.limit = qr(dat$X)$rank,
        rank.step = 2,
        pct = 0.98,
        lambda.multi.factor = 20,
        lambda.init = NULL,
        n.lambda = 20, 
        early.stopping = 1
      ),
      trace = F,
      quiet = T,
      track = F,
      step3 = T,
      use_warmstart = T,
      seed = NULL
    )
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = unsvd(fit.)
    fit.$beta = unsvd(fit.$beta)
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC-2")
    results$lambda.beta = fiti$hparams$lambda.beta
    results$lambda.M = fiti$hparams$lambda.M
    
    results <- c(results,
                          prepare_output(start_time, fit.$estimates, dat$O, dat$W, dat$beta, fit.$beta, dat$M, fit.$M, LogLik_SI))  
    
    results
  }
#----------------------------------------------------
CASMC_3a_Sim_Wrapper <-
  function(dat,
           max_cores = 20,
           LogLik_SI = NULL,
           ...) {
    start_time = Sys.time()
    learning_rate = 1 / sqrt(sum((t(dat$X) %*% dat$X)^2))
    fiti <- CASMC3_cv_beta(
      y_train = dat$fit_data$train,
      X = dat$X,
      y_valid = dat$fit_data$valid,
      W_valid = dat$fit_data$W_valid,
      y = dat$fit_data$Y,
      trace = 0,
      print.best = T,
      warm = NULL,
      quiet = F, learning.rate = learning_rate,
      early.stopping = 1,
      lambda.beta.grid = seq(0,10,length.out=20),
      max_cores = max_cores
    ) 
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC-3a")
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
