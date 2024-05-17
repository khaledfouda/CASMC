Mao_Sim_Wrapper <-
  function(dat,
           lambda.1_grid = seq(0, 1, length = 20),
           lambda.2_grid= seq(0.9, 0.1, length = 20),
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
      W = dat$W,
      n_folds = n_folds,
      lambda.1_grid = lambda.1_grid,
      lambda.2_grid = lambda.2_grid,
      alpha_grid = alpha_grid,#seq(0.992, 1, length = 5),
      seed = 2023,
      numCores = ncores,
      n1n2_optimized = TRUE,
      test_error = error_metric$rmse,
      theta_estimator = weight_function,
      sequential = FALSE
    )
    
    fit. <- fiti$fit
    results = list(model = "Mao")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$lambda.1 = fiti$best_parameters$lambda.1
    results$lambda.2 = fiti$best_parameters$lambda.2
    results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W ==
                                                                                   0])
    results$error.all = test_error(fit.$estimates, dat$O)
    results$error.M = test_error(fit.$M, dat$M)
    results$error.beta = test_error(fit.$beta, dat$beta)
    results$rank_M = qr(fit.$M)$rank
    results$rank_beta = qr(fit.$beta)$rank
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
    maxit = 300,
    seed= NULL
  )
  results = list(model = "SoftImpute")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$lambda.1 = NA
  results$lambda.2 = fit.$lambda
  results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W ==
                                                                              0])
  results$error.all = test_error(fit.$estimates, dat$O)
  results$error.M = NA
  results$error.beta = NA
  results$rank_M = fit.$rank_M
  results$rank_beta = NA
  results
}


CASMC_rank_Sim_Wrapper <-
  function(dat,
           max_cores=20,
           ...) {
    start_time = Sys.time()
    rank_x <- qr(dat$X)$rank
    fiti <- CASMC_cv_rank(
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
      maxit = 100,
      trace = T,
      print.best = F,
      quiet = FALSE,
      warm = NULL,
      rank_x = rank_x,
      r_min = 0,
      r_max = rank_x,
      track_r = F,
      max_cores = max_cores,
      seed = NULL
    )

    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = "CASMC_Rank_Restriction")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$lambda.1 = NA
    results$lambda.2 = fit.$lambda
    results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W == 0])
    results$error.all = test_error(fit.$estimates, dat$O)
    results$error.M = test_error(fit.$M, dat$M)
    results$error.beta = test_error(fit.$beta, dat$beta)
    results$rank_M = length(fit.$d)
    results$rank_beta = fiti$r
    results
  }




CASMC_L2_Sim_Wrapper <-
  function(dat,
           max_cores=20,
           maxit = 200,
           ...) {
    start_time = Sys.time()
    
    fiti <- CASMC_cv_L2(
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
      maxit = maxit,
      trace = FALSE,
      print.best = F,
      quiet = FALSE,
      warm = NULL,
      lambda.beta.grid = "default",
      track_beta = F,
      max_cores = max_cores,
      seed = NULL
    )
    
    fit. = fiti$fit
    # get estimates and validate
    fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
    fit.$estimates = fit.$M + dat$X %*% fit.$beta
    
    results = list(model = ifelse(maxit==2,"CASMC_L2_Beta_single_iter","CASMC_L2_Beta"))
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$lambda.1 = fiti$lambda.beta
    results$lambda.2 = fit.$lambda
    results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W == 0])
    results$error.all = test_error(fit.$estimates, dat$O)
    results$error.M = test_error(fit.$M, dat$M)
    results$error.beta = test_error(fit.$beta, dat$beta)
    results$rank_M = length(fit.$d)
    results$rank_beta = qr(fit.$beta)$rank
    results
  }






Naive_Sim_Wrapper <- function(dat, ...) {
  start_time = Sys.time()
  fiti <- naive_fit(dat$Y, dat$X) 
  results = list(model = "Naive")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$lambda.1 = NA
  results$lambda.2 = NA
  results$error.test = test_error(fiti$estimates[dat$W == 0], dat$O[dat$W == 0])
  results$error.all = test_error(fiti$estimates, dat$O)
  results$error.M = NA
  results$error.beta = NA 
  results$rank_M = qr(fiti$estimates)$rank
  results$rank_beta = NA
  results
}

