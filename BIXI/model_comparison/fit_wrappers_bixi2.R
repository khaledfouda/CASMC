


prepare_output_bixi <-
  function(start_time,
           X,
           estim.test,
           estim.train,
           obs.test,
           obs.train,
           beta.estim = NA,
           M.estim = NA,
           LogLik_SI = NA,
           test_error = error_metric$rmse) {
    list(
      time = round(as.numeric(
        difftime(Sys.time(), start_time, units = "secs")
      )),
      error.test = test_error(estim.test, obs.test),
      corr.test = cor(estim.test, obs.test),
      error.train = test_error(estim.train, obs.train),
      rank_M = tryCatch(
        qr(M.estim)$rank,
        error = function(x)
          NA
      ),
      rank_beta = tryCatch(
        qr(beta.estim)$rank,
        error = function(x)
          NA
      ),
      sparsity = tryCatch(
        sum(beta.estim == 0) /
          length(beta.estim),
        error = function(x)
          NA
      )
    ) -> results
    
    
    if (is.null(LogLik_SI)) {
      results$likelihood_ratio_index <- NA
      results$Cox_Snell_R2 <- NA
    } else{
      residuals <- obs.test - estim.test
      LogLik <- logLikelihood(residuals)
      n <- length(residuals)
      results$likelihood_ratio_index <-
        Likelihood_ratio_index(LogLik, LogLik_SI)
      results$Cox_Snell_R2 <- Cox_Snell_R2(LogLik, LogLik_SI, n)
    }
    
    
    tryCatch(
      apply(beta.estim, 1, summary) |> as.data.frame() |>
      t() |>
      as.data.frame() |>
      mutate(prop_non_zero = apply(beta.estim, 1, function(x)
        sum(x != 0) / length(x))) |>
      `rownames<-` (colnames(X)),
      error = function(x) NA) ->
      results$cov_summaries
    
    return(results)
  }

Mao_Bixi_Wrapper <-
  function(dat,
           lambda.1_grid = seq(0, 1, length = 20),
           lambda.2_grid = seq(0.9, 0.1, length = 20),
           alpha_grid = c(1),
           ncores = 1,
           # keep it > 1
           n_folds = 5,
           weight_function = Mao_weights$uniform,
           ...) {
    start_time = Sys.time()
    fiti <- Mao.cv(
      Y = dat$Y,
      X = dat$X,
      W = dat$masks$tr_val,
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
    results <- c(
      results,
      prepare_output_bixi(
        start_time,
        dat$X,
        fit.$estimates[dat$masks$test == 0],
        fit.$estimates[dat$masks$tr_val != 0],
        dat$splits$test@x,
        dat$splits$Y@x,
        fit.$beta,
        fit.$M
      )
    )
    results
  }

SImpute_Bixi_Wrapper <- function(dat, ...) {
  start_time = Sys.time()
  fit. <- simpute.cv(
    Y_train = dat$Y,
    y_valid = dat$splits$valid@x,
    W_valid = dat$masks$valid,
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
  results <- c(
    results,
    prepare_output_bixi(
      start_time,
      NULL,
      fit.$estimates[dat$masks$test == 0],
      fit.$estimates[dat$masks$tr_val != 0],
      dat$splits$test@x,
      dat$splits$Y@x,
      M.estim = fit.$estimates
    )
  )
  LogLik <-
    logLikelihood(dat$splits$test@x - fit.$estimates[dat$masks$test == 0])
  return(list(results = results, LogLik = LogLik))
}
#--------------------------------------------------------------------------------------
CASMC_0_Bixi_Wrapper <-
  function(dat,
           max_cores = 20,
           LogLik_SI = NULL,
           return_fit = FALSE,
           train_on_all = FALSE,
           ...) {
    start_time = Sys.time()
    Y_all <- NULL
    if(train_on_all) Y_all <- dat$splits$Y 
    
    fiti <- CASMC0_cv(
      y_train = dat$splits$train,
      X = dat$X,
      y_valid = dat$splits$valid@x,
      W_valid = dat$masks$valid,
      y = Y_all,
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
    fit.$Xbeta = dat$X %*% fit.$beta
    fit.$estimates = fit.$M + fit.$Xbeta
    
    print("hi")
    results = list(model = "CASMC-0")
    results$lambda.beta = fiti$lambda.beta
    results$lambda.M = fit.$lambda
    results <- c(
      results,
      prepare_output_bixi(
        start_time,
        dat$X,
        fit.$estimates[dat$masks$test == 0],
        fit.$estimates[dat$masks$tr_val != 0],
        dat$splits$test@x,
        dat$splits$Y@x,
        fit.$beta,
        fit.$M,
        LogLik_SI
      )
    )
    if(return_fit) return(list(fit=fit., results=results))
    results
  }
#-------


CASMC_2_Bixi_Wrapper <-
  function(dat,
           LogLik_SI = NULL,
           return_fit = FALSE,
           train_on_all = FALSE,
           ...) {
    start_time = Sys.time()
    Y_all <- NULL
    if(train_on_all) Y_all <- dat$splits$Y
    
    fiti <- CASMC2_cv2(
      y_train = dat$splits$train,
      X = dat$X,
      y_valid = dat$splits$valid@x,
      W_valid = dat$masks$valid,
      y = Y_all,
      error_function = error_metric$rmse,
      warm = NULL,
      M_cv_param = list(
        rank.init = 2,
        rank.limit = 30,
        rank.step = 2,
        pct = 0.98,
        lambda.factor = 1 / 4,
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
      quiet = T,
      trace = F,
      track = F,
      step3 = T,
      use_warmstart = TRUE,
      seed = NULL,
    )
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = unsvd(fit.)
    fit.$beta = unsvd(fit.$beta)
    fit.$Xbeta = dat$X %*% fit.$beta
    fit.$estimates = fit.$M + fit.$Xbeta
    
    results = list(model = "CASMC-2")
    results$lambda.beta = fiti$hparams$lambda.beta
    results$lambda.M = fiti$hparams$lambda.M
    results <- c(
      results,
      prepare_output_bixi(
        start_time,
        dat$X,
        fit.$estimates[dat$masks$test == 0],
        fit.$estimates[dat$masks$tr_val != 0],
        dat$splits$test@x,
        dat$splits$Y@x,
        fit.$beta,
        fit.$M,
        LogLik_SI
      )
    )
    if(return_fit) return(list(fit=fit., results=results))
    results
  }
#----------------------------------------------------
CASMC_3a_Bixi_Wrapper <-
  function(dat,
           max_cores = 20,
           LogLik_SI = NULL,
           return_fit = FALSE,
           train_on_all = FALSE,
           ...) {
    start_time = Sys.time()
    learning_rate = 1 / sqrt(sum((t(dat$X) %*% dat$X) ^ 2))
    Y_all <- NULL
    if(train_on_all) Y_all <- dat$splits$Y
    
    fiti <- CASMC3_cv_beta(
      y_train = dat$splits$train,
      X = dat$X,
      y_valid = dat$splits$valid@x,
      W_valid = dat$masks$valid,
      y = Y_all,
      trace = 0,
      print.best = T,
      warm = NULL,
      quiet = F,
      learning.rate = learning_rate,
      early.stopping = 1,
      lambda.beta.grid = seq(0, 10, length.out = 20),
      max_cores = max_cores
    )
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$Xbeta = dat$X %*% fit.$beta
    fit.$estimates = fit.$M + fit.$Xbeta
    
    results = list(model = "CASMC-3a")
    results$lambda.beta = fiti$hparams$lambda.beta
    results$lambda.M = fiti$hparams$lambda.M
    results <- c(
      results,
      prepare_output_bixi(
        start_time,
        dat$X,
        fit.$estimates[dat$masks$test == 0],
        fit.$estimates[dat$masks$tr_val != 0],
        dat$splits$test@x,
        dat$splits$Y@x,
        fit.$beta,
        fit.$M,
        LogLik_SI
      )
    )
    if(return_fit) return(list(fit=fit., results=results))
    results
  }
#------
#----------------------------------------------
Naive_Bixi_Wrapper <- function(dat, ...) {
  start_time = Sys.time()
  fit. <- naive_fit(dat$Y, dat$X)
  results = list(model = "Naive")
  results$lambda.beta = NA
  results$lambda.M = NA
  results <- c(
    results,
    prepare_output_bixi(
      start_time,
      dat$X,
      fit.$estimates[dat$masks$test == 0],
      fit.$estimates[dat$masks$tr_val != 0],
      dat$splits$test@x,
      dat$splits$Y@x,
      fit.$beta,
      fit.$M
    )
  )
  
  results
}
