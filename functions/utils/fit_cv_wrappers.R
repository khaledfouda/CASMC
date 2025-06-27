#------------------------------------------------------------------------------#
# Prepare a standardized output list for model fits
#------------------------------------------------------------------------------#
prepare_output <- function(
    start_time,
    estimates,
    obs,
    mask,
    beta        = NA,
    beta.estim  = NA,
    M           = NA,
    M.estim     = NA,
    LogLik_SI   = NA,
    test_error  = utils$error_metric$rmse
) {
  # Split train vs test
  estim.test  <- estimates[mask == 0]
  estim.train <- estimates[mask != 0]
  obs.test    <- obs[mask == 0]
  obs.train   <- obs[mask != 0]
  
  # Core metrics
  results <- list(
    time       = round(as.numeric(difftime(
      Sys.time(), start_time, units = "secs"
    ))),
    error.test = test_error(estim.test, obs.test),
    corr.test  = cor(estim.test, obs.test),
    error.train = test_error(estim.train, obs.train),
    error.M    = tryCatch(
      test_error(M.estim, M),
      error = function(e) NA
    ),
    error.beta = tryCatch(
      test_error(beta.estim, beta),
      error = function(e) NA
    ),
    rank_M     = tryCatch(
      qr(M.estim)$rank,
      error = function(e) NA
    ),
    rank_beta  = tryCatch(
      qr(beta.estim)$rank,
      error = function(e) NA
    ),
    sparse_in_sparse = tryCatch(
      sum(beta == 0 & beta.estim == 0) /
        (sum(beta == 0) + 1e-17),
      error = function(e) NA
    ),
    nonsparse_in_nonsparse = tryCatch(
      sum(beta != 0 & beta.estim != 0) /
        (sum(beta != 0) + 1e-17),
      error = function(e) NA
    ),
    sparse_all = tryCatch(
      sum(beta.estim == 0) / length(beta.estim),
      error = function(e) NA
    )
  )
  # 
  # # Add pseudo-R² only if SI log-lik was provided
  # if (is.null(LogLik_SI)) {
  #   results$likelihood_ratio_index <- NA
  #   results$Cox_Snell_R2          <- NA
  # } else {
  #   residuals <- obs.test - estim.test
  #   LogLik    <- utils$logLikelihood(residuals)
  #   n         <- length(residuals)
  #   
  #   results$likelihood_ratio_index <- utils$Likelihood_ratio_index(
  #     LogLik, LogLik_SI
  #   )
  #   results$Cox_Snell_R2 <- utils$Cox_Snell_R2(
  #     LogLik, LogLik_SI, n
  #   )
  # }
  
  results
}

#------------------------------------------------------------------------------#
# Wrapper for the “Mao” method with cross-validation
#------------------------------------------------------------------------------#
Mao_Sim_Wrapper <- function(
    dat,
    lambda_1_grid    = seq(0,   1,   length = 20),
    lambda_2_grid    = seq(0.9, 0.1, length = 20),
    alpha_grid       = 1,
    ncores           = 1,
    n_folds          = 5,
    weight_function  = Mao_weights$uniform,
    LogLik_SI        = NULL,
    ...
) {
  start_time <- Sys.time()
  
  cv_out <- Mao.cv(
    Y               = dat$Y,
    X               = dat$fit_data$Xq,
    W               = dat$W,
    n_folds         = n_folds,
    lambda_1_grid   = lambda_1_grid,
    lambda_2_grid   = lambda_2_grid,
    alpha_grid      = alpha_grid,
    seed            = 2023,
    numCores        = ncores,
    n1n2_optimized  = TRUE,
    test_error      = utils$error_metric$rmse,
    theta_estimator = weight_function,
    sequential      = FALSE
  )
  
  fit <- cv_out$fit
  results <- list(model = "Mao")
  results$lambda_beta <- cv_out$best_parameters$lambda_1
  results$lambda_M    <- cv_out$best_parameters$lambda_2
  
  results <- c(
    results,
    prepare_output(
      start_time = start_time,
      estimates  = fit$estimates,
      obs        = dat$O,
      mask       = dat$W,
      beta       = dat$beta,
      beta.estim = fit$beta,
      M          = dat$M,
      M.estim    = fit$M,
      LogLik_SI  = LogLik_SI
    )
  )
  
  results
}

#------------------------------------------------------------------------------#
# Wrapper for SoftImpute via simpute.cv
#------------------------------------------------------------------------------#
SImpute_Sim_Wrapper <- function(dat, ...) {
  start_time <- Sys.time()
  
  fit <- simpute.cv(
    Y_train   = as.matrix(dat$fit_data$train),
    y_valid   = dat$fit_data$valid,
    W_valid   = dat$fit_data$W_valid,
    y         = dat$Y,
    n.lambda  = 20,
    trace     = FALSE,
    print.best= FALSE,
    tol       = 5,
    thresh    = 1e-6,
    rank.init = 2,
    rank.limit= 30,
    rank.step = 2,
    maxit     = 600,
    seed      = NULL
  )
  
  results <- list(model = "SoftImpute")
  results$lambda_beta <- NA
  results$lambda_M    <- fit$lambda
  
  results <- c(
    results,
    prepare_output(
      start_time = start_time,
      estimates  = fit$estimates,
      obs        = dat$O,
      mask       = dat$W,
      M.estim    = fit$estimates
    )
  )
  
  # Return both results and log-lik for further comparison
  LogLik <- utils$logLikelihood(
    dat$O[dat$W == 0] - fit$estimates[dat$W == 0]
  )
  
  list(results = results, LogLik = LogLik)
}

#------------------------------------------------------------------------------#
# Wrapper for CASMC with Lasso penalty
#------------------------------------------------------------------------------#
CAMC_Sim_Wrapper <- function(
    dat,
    max_cores   = 20,
    LogLik_SI   = NULL,
    hpar        = CAMC_Lasso_hparams,
    verbose       = 1,
    return_fit  = FALSE,
    ...
) {
  start_time <- Sys.time()
  
  cv_out <- CAMC_Lasso_cv(
    y_train    = dat$fit_data$train,
    X          = dat$fit_data$Xq,
    y_valid    = dat$fit_data$valid,
    W_valid    = dat$fit_data$W_valid,
    y          = dat$fit_data$Y,
    hpar       = hpar,
    verbose    = verbose,
    max_cores  = max_cores
  )
  
  fit <- cv_out$fit
  fit$M         <- fit$u %*% (fit$d * t(fit$v))
  fit$estimates <- fit$M + dat$fit_data$Xq %*% fit$beta
  
  results <- list(model = "CASMC-Lasso")
  results$lambda_beta <- cv_out$hparams$lambda_beta
  results$lambda_M    <- cv_out$hparams$lambda_M
  
  results <- c(
    results,
    prepare_output(
      start_time = start_time,
      estimates  = fit$estimates,
      obs        = dat$O,
      mask       = dat$W,
      beta       = dat$beta,
      beta.estim = fit$beta,
      M          = dat$M,
      M.estim    = fit$M,
      LogLik_SI  = LogLik_SI
    )
  )
  
  if (return_fit) {
    return(list(results = results, fit = cv_out))
  }
  
  results
}

#------------------------------------------------------------------------------#
# Naive baseline wrapper
#------------------------------------------------------------------------------#
Naive_Sim_Wrapper <- function(dat, ...) {
  start_time <- Sys.time()
  fit        <- naive_fit(dat$Y, dat$fit_data$Xq)
  
  results <- list(model = "Naive")
  results$lambda_beta <- NA
  results$lambda_M    <- NA
  
  results <- c(
    results,
    prepare_output(
      start_time = start_time,
      estimates  = fit$estimates,
      obs        = dat$O,
      mask       = dat$W,
      beta       = dat$beta,
      beta.estim = fit$beta,
      M          = dat$M,
      M.estim    = fit$M
    )
  )
  
  results
}
