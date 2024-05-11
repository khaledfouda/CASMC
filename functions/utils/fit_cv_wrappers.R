compare_Mao <-
  function(gen.dat,
           lambda.1_grid,
           lambda.2_grid,
           alpha_grid,
           ncores = 2,
           # keep it > 1
           n_folds = 5,
           weight_function = MaoBinomalWeights,
           ...) {
    start_time = Sys.time()
    cv.out <- Mao.cv(
      gen.dat$O,
      gen.dat$X,
      gen.dat$Y,
      gen.dat$W,
      n_folds = n_folds,
      lambda.1_grid = lambda.1_grid,
      lambda.2_grid = lambda.2_grid,
      alpha_grid = alpha_grid,
      numCores = ncores,
      n1n2_optimized = TRUE,
      #keep it at TRUE
      theta_estimator = weight_function
    )
    mao.out <-
      Mao.fit(
        gen.dat$Y,
        gen.dat$X,
        gen.dat$W,
        cv.out$best_parameters$lambda.1,
        cv.out$best_parameters$lambda.2,
        cv.out$best_parameters$alpha,
        theta_estimator = weight_function,
        n1n2_optimized = TRUE
      )
    
    results = list(model = "Mao")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$alpha = cv.out$best_parameters$alpha
    results$lambda.1 = cv.out$best_parameters$lambda.1
    results$lambda.2 = cv.out$best_parameters$lambda.2
    results$error.test = test_error(mao.out$estimates[gen.dat$W == 0], gen.dat$O[gen.dat$W ==
                                                                                   0])
    results$error.all = test_error(mao.out$estimates, gen.dat$O)
    results$error.M = test_error(mao.out$M, gen.dat$M)
    results$error.beta = test_error(mao.out$beta, gen.dat$beta)
    results$rank = mao.out$rank
    results
  }

compare_softImpute <- function(gen.dat, valid.dat, ...) {
  start_time = Sys.time()
  sout <- simpute.cv(
    valid.dat$Y_train,
    gen.dat$Y,
    trace = FALSE,
    rank.limit = 30,
    print.best = FALSE,
    rank.step = 4
  )
  results = list(model = "SoftImpute")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$alpha = NA
  results$lambda.1 = NA
  results$lambda.2 = sout$lambda
  results$error.test = test_error(sout$estimates[gen.dat$W == 0], gen.dat$O[gen.dat$W ==
                                                                              0])
  results$error.all = test_error(sout$estimates, gen.dat$O)
  results$error.M = NA
  results$error.beta = NA
  results$rank = sout$rank_M
  results
}

compare_CAMC_holdout <-
  function(gen.dat,
           valid.dat,
           lambda.1_grid,
           rank.step,
           rank.limit,
           n.lambda,
           ...) {
    start_time = Sys.time()
    
    fiti <-
      CAMC_cv_holdout(
        valid.dat$Y_train,
        gen.dat$X,
        valid.dat$W_valid,
        valid.dat$Y_valid,
        trace = FALSE,
        rank.limit = rank.limit,
        print.best = FALSE,
        rank.step  = rank.step,
        n.lambda = n.lambda,
        type = "als",
        quiet = TRUE,
        tol = 2,
        lambda.1_grid = lambda.1_grid
      )
    sout <- fiti$fit2
    
    
    results = list(model = "CAMC_holdout")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$alpha = NA
    results$lambda.1 = sout$lambda1
    results$lambda.2 = sout$lambda2
    results$error.test = test_error(sout$estimates[gen.dat$W == 0], gen.dat$O[gen.dat$W ==
                                                                                0])
    results$error.all = test_error(sout$estimates, gen.dat$O)
    results$error.M = test_error(sout$M, gen.dat$M)
    results$error.beta = test_error(sout$beta, gen.dat$beta)
    results$rank = sout$rank_O
    results
  }

compare_CAMC_kfold <-
  function(gen.dat,
           lambda.1_grid,
           n_folds,
           rank.step,
           rank.limit,
           n.lambda,
           ...) {
    start_time = Sys.time()
    fiti <- CAMC_cv_kfold(
      gen.dat$Y,
      gen.dat$X,
      gen.dat$W,
      n_folds = n_folds,
      trace = FALSE,
      rank.limit = rank.limit,
      print.best = FALSE,
      lambda.1_grid = lambda.1_grid,
      rank.step = rank.step,
      n.lambda = n.lambda,
      type = "als",
      tol = 2
    )
    sout <- fiti$fit2
    results = list(model = "CAMC_kfold")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$alpha = NA
    results$lambda.1 = sout$lambda1
    results$lambda.2 = sout$lambda2
    results$error.test = test_error(sout$estimates[gen.dat$W == 0], gen.dat$O[gen.dat$W ==
                                                                                0])
    results$error.all = test_error(sout$estimates, gen.dat$O)
    results$error.M = test_error(sout$M, gen.dat$M)
    results$error.beta = test_error(sout$beta, gen.dat$beta)
    results$rank = sout$rank_O
    results
  }




compare_CASMC_holdout <-
  function(gen.dat,
           valid.dat,
           splr.dat,
           rank.step,
           rank.limit,
           n.lambda,
           max_cores=20,
           ...) {
    start_time = Sys.time()
    
    y_train = valid.dat$Y_train
    y_train[y_train == 0] = NA
    y_train = as(y_train, "Incomplete")
    
    y = gen.dat$Y
    y[y == 0] = NA
    y = as(y, "Incomplete")
    
    best_fit = CASMC_cv_holdout_with_reg(
      y_train,
      splr.dat,
      valid.dat$Y_valid,
      valid.dat$W_valid,
      #r_min = 0,
      y = y,
      trace = F,
      max_cores = max_cores,
      thresh = 1e-6,
      n.lambda = n.lambda,
      rank.limit = rank.limit,
      maxit = 200,
      rank.step = rank.step,
      print.best = FALSE
    )
    print(best_fit$r)
    
    fit1 = best_fit$fit
    #fit2 = best_fit$fit2
    sout = best_fit
    # get estimates and validate
    sout$M = fit1$u %*% (fit1$d * t(fit1$v))
    sout$beta =  t(fit1$beta)#fit1$Beta$u %*% (fit1$Beta$d * t(fit1$Beta$v))
    sout$estimates = sout$M + splr.dat$X %*% t(sout$beta)
    
    results = list(model = "CASMC_holdout")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$alpha = NA
    results$lambda.1 = sout$lambda.beta
    results$lambda.2 = sout$lambda
    results$error.test = test_error(sout$estimates[gen.dat$W == 0], gen.dat$O[gen.dat$W ==
                                                                                0])
    results$error.all = test_error(sout$estimates, gen.dat$O)
    results$error.M = test_error(sout$M, gen.dat$M)
    results$error.beta = test_error(t(sout$beta), gen.dat$beta)
    results$rank = qr(sout$estimates)$rank
    results
  }

compare_CASMC_kfold <-
  function(gen.dat,
           splr.dat,
           n_folds,
           rank.step,
           rank.limit,
           n.lambda,
           ...) {
    start_time = Sys.time()
    
    best_fit = CASMC_cv_kfold_v2(
      gen.dat$Y,
      splr.dat,
      gen.dat$W,
      trace = F,
      print.best = F,
      rank.limit = rank.limit,
      n.lambda = n.lambda,
      maxit = 200,
      n_folds = n_folds,
      rank.step = rank.step
    )
    
    
    fit1 = best_fit$fit
    sout = best_fit
    # get estimates and validate
    sout$M = fit1$u %*% (fit1$d * t(fit1$v))
    sout$beta =  fit1$beta
    sout$estimates = sout$M + splr.dat$X %*% t(sout$beta)
    
    
    results = list(model = "CASMC_kfold")
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results$alpha = NA
    results$lambda.1 = NA
    results$lambda.2 = sout$lambda
    results$error.test = test_error(sout$estimates[gen.dat$W == 0], gen.dat$O[gen.dat$W ==
                                                                                0])
    results$error.all = test_error(sout$estimates, gen.dat$O)
    results$error.M = test_error(sout$M, gen.dat$M)
    results$error.beta = test_error(t(sout$beta), gen.dat$beta)
    results$rank = qr(sout$estimates)$rank
    results
  }

compare_naive <- function(gen.dat, ...) {
  start_time = Sys.time()
  estimates = naive_MC(gen.dat$Y)
  results = list(model = "Naive")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$alpha = NA
  results$lambda.1 = NA
  results$lambda.2 = NA
  results$error.test = test_error(estimates[gen.dat$W == 0], gen.dat$O[gen.dat$W == 0])
  results$error.all = test_error(estimates, gen.dat$O)
  results$error.M = NA
  results$error.beta = NA
  results$rank = qr(estimates)$rank
  results
}


compare_cmf <- function(gen.dat, splr.dat, ...) {
  start_time = Sys.time()
  
  library(cmfrec)
  library(Matrix)
  library(MatrixExtra)
  
  X_test = gen.dat$O
  X_test[gen.dat$W == 1] = NA
  X_test = as(X_test, "Incomplete")
  X_test <- as.coo.matrix(X_test)
  
  X_train <- gen.dat$Y
  X_train[gen.dat$W == 0] = NA
  X_train <- as(X_train, "Incomplete")
  X_train <- as.coo.matrix(X_train)
  
  model.w.sideinfo <- CMF(
    X_train,
    U = splr.dat$X,
    k = 25,
    lambda = 0.1,
    scale_lam = TRUE,
    niter = 30,
    use_cg = FALSE,
    include_all_X = FALSE,
    w_main = 0.75,
    w_user = 0.5,
    w_item = 0.5,
    w_implicit = 0.5,
    center_U = FALSE,
    center_I = FALSE,
    nthreads = 6,
    verbose = FALSE
  )
  pred_side_info <- predict(model.w.sideinfo, X_test)
  #print_rmse(X_test, pred_side_info, "model with side info")
  results = list(model = "CMF")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  detach("package:cmfrec", unload = TRUE)
  detach("package:MatrixExtra", unload = TRUE)
  source("./code_files/import_lib.R")
  
  results$alpha = NA
  results$lambda.1 = NA
  results$lambda.2 = NA
  results$error.test = test_error(pred_side_info@x, X_test@x)
  results$error.all = NA #test_error(sout$estimates, gen.dat$O)
  results$error.M = NA
  results$error.beta = NA
  results$rank = NA#sout$rank_M
  results
}