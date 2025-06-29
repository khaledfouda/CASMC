#------------------------------------------------------------------------------#
# Prepare standardized output for BIXI-style wrappers
#------------------------------------------------------------------------------#
prepare_output_bixi <- function(
    start_time,
    X,
    estim.test,
    estim.train,
    obs.test,
    obs.train,
    beta.estim  = NA,
    M.estim     = NA,
    #LogLik_SI   = NA,
    test_error  = utils$error_metric$rmse
) {
  # Core metrics
  results <- list(
    time        = round(as.numeric(
      difftime(Sys.time(), start_time, units = "secs")
    )),
    error.test  = test_error(estim.test, obs.test),
    corr.test   = cor(estim.test, obs.test),
    error.train = test_error(estim.train, obs.train),
    rank_M      = tryCatch(
      qr(M.estim)$rank,
      error = function(e) NA
    ),
    rank_beta   = tryCatch(
      qr(beta.estim)$rank,
      error = function(e) NA
    ),
    sparsity    = tryCatch(
      sum(beta.estim == 0) / length(beta.estim),
      error = function(e) NA
    )
  )
  
  # # Pseudo-R² if baseline log-likelihood provided
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
  
  # Covariate coefficient summaries
  results$cov_summaries <- tryCatch({
    apply(beta.estim, 1, summary) |>
      as.data.frame() |>
      t() |>
      as.data.frame() |>
      dplyr::mutate(
        prop_non_zero = apply(beta.estim, 1, function(x)
          sum(x != 0) / length(x)
        )
      ) |>
      `rownames<-`(colnames(X))
  }, error = function(e) NA)
  
  results
}

#------------------------------------------------------------------------------#
# Wrapper for Mao method adapted to BIXI data
#------------------------------------------------------------------------------#
Mao_Bixi_Wrapper <- function(
    dat,
    lambda_1_grid   = seq(0,   1,   length = 20),
    lambda_2_grid   = seq(0.9, 0.1, length = 20),
    alpha_grid      = 1,
    ncores          = 1,
    n_folds         = 5,
    weight_function = Mao_weights$uniform,
    ...
) {
  start_time <- Sys.time()
  
  cv_out <- Mao.cv(
    Y               = dat$Y,
    X               = dat$X,
    W               = dat$masks$tr_val,
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
    prepare_output_bixi(
      start_time  = start_time,
      X           = dat$X,
      estim.test  = fit$estimates[dat$masks$test == 0],
      estim.train = fit$estimates[dat$masks$tr_val != 0],
      obs.test    = dat$splits$test@x,
      obs.train   = dat$splits$Y@x,
      beta.estim  = fit$beta,
      M.estim     = fit$M
    )
  )
  
  results
}

#------------------------------------------------------------------------------#
# Wrapper for SoftImpute (simpute.cv) adapted to BIXI data
#------------------------------------------------------------------------------#
SImpute_Bixi_Wrapper <- function(dat, ...) {
  start_time <- Sys.time()
  
  fit <- simpute.cv(
    Y_train   = dat$Y,
    y_valid   = dat$splits$valid@x,
    W_valid   = dat$masks$valid,
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
    prepare_output_bixi(
      start_time  = start_time,
      X           = NULL,
      estim.test  = fit$estimates[dat$masks$test == 0],
      estim.train = fit$estimates[dat$masks$tr_val != 0],
      obs.test    = dat$splits$test@x,
      obs.train   = dat$splits$Y@x,
      M.estim     = fit$estimates
    )
  )
  
  LogLik <- utils$logLikelihood(
    dat$splits$test@x - fit$estimates[dat$masks$test == 0]
  )
  list(results = results, LogLik = LogLik)
}

#------------------------------------------------------------------------------#
# Wrapper for CAMC with Ridge penalty adapted to BIXI data
#------------------------------------------------------------------------------#
CAMC_0_Bixi_Wrapper <- function(
    dat,
    max_cores    = 20,
    LogLik_SI    = NULL,
    return_fit   = FALSE,
    train_on_all = FALSE,
    hpar         = CAMC_Ridge_hparams,
    ...
) {
  start_time <- Sys.time()
  Y_all      <- if (train_on_all) dat$splits$Y else NULL
  
  cv_out <- CAMC_Ridge_cv(
    y_train    = dat$splits$train,
    X          = dat$X,
    y_valid    = dat$splits$valid@x,
    W_valid    = dat$masks$valid,
    y          = Y_all,
    hpar       = hpar,
    thresh     = 1e-6,
    maxit      = 300,
    trace      = TRUE,
    track      = TRUE,
    print.best = FALSE,
    quiet      = FALSE,
    warm       = NULL,
    max_cores  = max_cores,
    seed       = NULL
  )
  
  fit <- cv_out$fit
  fit$M         <- fit$u %*% (fit$d * t(fit$v))
  fit$Xbeta     <- dat$X %*% fit$beta
  fit$estimates <- fit$M + fit$Xbeta
  
  results <- list(model = "CAMC-0")
  results$lambda_beta <- cv_out$lambda_beta
  results$lambda_M    <- fit$lambda
  
  results <- c(
    results,
    prepare_output_bixi(
      start_time  = start_time,
      X           = dat$X,
      estim.test  = fit$estimates[dat$masks$test == 0],
      estim.train = fit$estimates[dat$masks$tr_val != 0],
      obs.test    = dat$splits$test@x,
      obs.train   = dat$splits$Y@x,
      beta.estim  = fit$beta,
      M.estim     = fit$M,
      LogLik_SI   = LogLik_SI
    )
  )
  
  if (return_fit) {
    return(list(fit = fit, results = results))
  }
  
  results
}

#------------------------------------------------------------------------------#
# Wrapper for CAMC with Nuclear-norm penalty adapted to BIXI data
#------------------------------------------------------------------------------#
CAMC_2_Bixi_Wrapper <- function(
    dat,
    LogLik_SI    = NULL,
    return_fit   = FALSE,
    train_on_all = FALSE,
    hpar         = CAMC_Nuclear_hparams,
    ...
) {
  start_time <- Sys.time()
  Y_all      <- if (train_on_all) dat$splits$Y else NULL
  
  cv_out <- CAMC_Nuclear_cv(
    y_train       = dat$splits$train,
    X             = dat$X,
    y_valid       = dat$splits$valid@x,
    W_valid       = dat$masks$valid,
    y             = Y_all,
    hpar          = hpar,
    warm          = NULL,
    quiet         = TRUE,
    trace         = FALSE,
    track         = FALSE,
    step3         = TRUE,
    use_warmstart = TRUE,
    seed          = NULL
  )
  
  fit <- cv_out$fit
  fit$M         <- utils$unsvd(fit)
  fit$beta      <- utils$unsvd(fit$beta)
  fit$Xbeta     <- dat$X %*% fit$beta
  fit$estimates <- fit$M + fit$Xbeta
  
  results <- list(model = "CAMC-2")
  results$lambda_beta <- cv_out$hparams$lambda_beta
  results$lambda_M    <- cv_out$hparams$lambda_M
  
  results <- c(
    results,
    prepare_output_bixi(
      start_time  = start_time,
      X           = dat$X,
      estim.test  = fit$estimates[dat$masks$test == 0],
      estim.train = fit$estimates[dat$masks$tr_val != 0],
      obs.test    = dat$splits$test@x,
      obs.train   = dat$splits$Y@x,
      beta.estim  = fit$beta,
      M.estim     = fit$M,
      LogLik_SI   = LogLik_SI
    )
  )
  
  if (return_fit) {
    return(list(fit = fit, results = results))
  }
  
  results
}

#------------------------------------------------------------------------------#
# Wrapper for CAMC with Lasso penalty adapted to BIXI data
#------------------------------------------------------------------------------#
CAMC_3a_Bixi_Wrapper <- function(
    dat,
    max_cores    = 20,
    LogLik_SI    = NULL,
    return_fit   = FALSE,
    train_on_all = FALSE,
    hpar         = CAMC_Lasso_hparams,
    ...
) {
  start_time <- Sys.time()
  Y_all      <- if (train_on_all) dat$splits$Y else NULL
  
  cv_out <- CAMC_Lasso_cv(
    y_train    = dat$splits$train,
    X          = dat$X,
    y_valid    = dat$splits$valid@x,
    W_valid    = dat$masks$valid,
    y          = Y_all,
    hpar       = hpar,
    trace      = 3,
    print.best = TRUE,
    warm       = NULL,
    quiet      = FALSE,
    max_cores  = max_cores
  )
  
  fit <- cv_out$fit
  fit$M         <- fit$u %*% (fit$d * t(fit$v))
  fit$Xbeta     <- dat$X %*% fit$beta
  fit$estimates <- fit$M + fit$Xbeta
  
  results <- list(model = "CAMC-3a")
  results$lambda_beta <- cv_out$hparams$lambda_beta
  results$lambda_M    <- cv_out$hparams$lambda_M
  
  results <- c(
    results,
    prepare_output_bixi(
      start_time  = start_time,
      X           = dat$X,
      estim.test  = fit$estimates[dat$masks$test == 0],
      estim.train = fit$estimates[dat$masks$tr_val != 0],
      obs.test    = dat$splits$test@x,
      obs.train   = dat$splits$Y@x,
      beta.estim  = fit$beta,
      M.estim     = fit$M,
      LogLik_SI   = LogLik_SI
    )
  )
  
  if (return_fit) {
    return(list(fit = fit, results = results))
  }
  
  results
}

#------------------------------------------------------------------------------#
# Naive baseline wrapper for BIXI data
#------------------------------------------------------------------------------#
Naive_Bixi_Wrapper <- function(dat, ...) {
  start_time <- Sys.time()
  fit        <- naive_fit(dat$Y, dat$X)
  
  results <- list(model = "Naive")
  results$lambda_beta <- NA
  results$lambda_M    <- NA
  
  results <- c(
    results,
    prepare_output_bixi(
      start_time  = start_time,
      X           = dat$X,
      estim.test  = fit$estimates[dat$masks$test == 0],
      estim.train = fit$estimates[dat$masks$tr_val != 0],
      obs.test    = dat$splits$test@x,
      obs.train   = dat$splits$Y@x,
      beta.estim  = fit$beta,
      M.estim     = fit$M
    )
  )
  
  results
}