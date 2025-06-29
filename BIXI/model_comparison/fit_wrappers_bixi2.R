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
    sequential      = FALSE,
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
    sequential      = sequential
  )
  
  fit <- cv_out$fit
  grid_size = ifelse(sequential == TRUE,
                     length(alpha_grid) + length(lambda_2_grid),
                     length(alpha_grid) * length(lambda_2_grid))
  grid_size = grid_size + length(lambda_1_grid)
  grid_size = grid_size * n_folds
  results <- list(model = "Mao", 
                  grid_size = grid_size)
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
SImpute_Bixi_Wrapper <- function(dat,
                                 hpar = CAMC_Lasso_hparams,
                                 ...) {
  start_time <- Sys.time()
  
  fit <- simpute.cv(
    Y_train   = dat$Y,
    y_valid   = dat$splits$valid@x,
    W_valid   = dat$masks$valid,
    y         = dat$Y,
    n.lambda  = hpar$M$n.lambda,
    trace     = FALSE,
    print.best= FALSE,
    tol       = 5,
    thresh    = 1e-6,
    rank.init = hpar$M$rank.init,
    rank.limit= hpar$M$rank.limit,
    rank.step = hpar$M$rank.step,
    maxit     = 600,
    seed      = NULL
  )
  
  grid_size <- paste0("M(", hpar$M$n.lambda, ")")
  results <- list(model = "SoftImpute",
                  grid_size = grid_size)
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
  
  # LogLik <- utils$logLikelihood(
  #   dat$splits$test@x - fit$estimates[dat$masks$test == 0]
  # )
  # list(results = results, LogLik = LogLik)
  return(results)
}



#------------------------------------------------------------------------------#
# Wrapper for CAMC with Lasso penalty adapted to BIXI data
#------------------------------------------------------------------------------#
CAMC_Bixi_Wrapper <- function(
    dat,
    max_cores    = 20,
    #LogLik_SI    = NULL,
    return_fit   = FALSE,
    train_on_all = FALSE,
    verbose      = 1,
    hpar         = CAMC_Lasso_hparams,
    ...
) {
  start_time <- Sys.time()
  Y_all      <- if (train_on_all) dat$splits$Y else NULL
  
  cv_out <- CAMC_Lasso_cv(
    y_train    = dat$splits$train,
    X          = dat$splits$Xq,
    y_valid    = dat$splits$valid@x,
    W_valid    = dat$masks$valid,
    y          = Y_all,
    hpar       = hpar,
    verbose      = verbose,
    max_cores  = max_cores
  )
  
  fit <- cv_out$fit
  fit$Rbeta     <- fit$beta
  fit$beta      <- MASS::ginv(dat$splits$Xr) %*% fit$Rbeta
  fit$beta      <- round(fit$beta, 10) # set small values to 0
  fit$M         <- fit$u %*% (fit$d * t(fit$v))
  fit$Xbeta     <- dat$splits$Xq %*% fit$Rbeta #dat$X %*% fit$beta
  fit$estimates <- fit$M + fit$Xbeta
  
  grid_size <- paste0("M(", hpar$M$n.lambda, ")*", hpar$beta$n.lambda)
  results <- list(model = "CAMC",
                  grid_size = grid_size)
  results$lambda_beta <- cv_out$lambda_beta
  results$lambda_M    <- cv_out$lambda_M
  
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
      #LogLik_SI   = LogLik_SI
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
  
  results <- list(model = "Naive", grid_size = 0)
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