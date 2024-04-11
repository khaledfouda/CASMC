print_rmse <- function(X_test, X_hat, model_name) {
  rmse <- sqrt(mean((X_test@x - X_hat@x) ^ 2))
  cat(sprintf("RMSE for %s is: %.4f\n", model_name, rmse))
}
#---------------
apply_to_sim_dat <- function(dat) {
  # if("package:cmfrec" %in% search()) detach("package:cmfrec", unload=TRUE)
  # if("package:MatrixExtra" %in% search()) detach("package:MatrixExtra", unload=TRUE)
  # if("package:Matrix" %in% search()) detach("package:Matrix", unload=TRUE,force = T)
  # library(Matrix)
  # source("./code_files/import_lib.R")
  
  dat$X_r = reduced_hat_decomp(dat$X, 1e-2)
  dat$train.inc = dat$Y
  dat$train.inc[dat$train.inc == 0] = NA
  dat$train.mat = dat$train.inc
  dat$train.inc <- as(dat$train.inc, "Incomplete")
  
  dat$test.inc = dat$O * (1 - dat$W)
  dat$test.inc[dat$test.inc == 0] = NA
  dat$test.inc <- as(dat$test.inc, "Incomplete")
  
  #------------------------------------------------
  
  
  W_valid <- matrix.split.train.test(dat$W, testp = 0.2)
  Y_train = (dat$Y * W_valid)
  Y_valid = dat$Y[W_valid == 0]
  
  
  start_time <- Sys.time()
  best_fit = CASMC_cv_holdout_v2(
    Y_train,
    dat$X_r,
    Y_valid,
    W_valid,
    dat$Y,
    trace = F,
    maxit = 100,
    thresh = 1e-6,
    n.lambda = 20,
    rank.step = 2,
    rank.limit = 30
  )
  fit1 = best_fit$fit
  beta =  fit1$beta
  M = fit1$u %*% (fit1$d * t(fit1$v))
  A = M + dat$X_r$X %*% t(beta)
  #A = revertBiScaledMatrix(as.matrix(A), dat$biScale)
  #qr(A)$rank
  preds = A[dat$W == 0]
  time_casmc <-
    as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  # RMSE_error(preds, dat$O[dat$W==0])
  #----------------------------------------------------------------------------------
  # soft-impute
  start_time = Sys.time()
  sout <- simpute.cv(
    Y_train,
    dat$Y,
    trace = FALSE,
    rank.limit = 30,
    print.best = FALSE,
    rank.step = 2
  )
  preds_simpute <- sout$estimates[dat$W == 0]
  
  time_softImpute = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  #---------------------------------------------------------------------------------
  # Naive
  start_time = Sys.time()
  estimates = naive_MC(dat$Y)
  preds_naive = estimates[dat$W == 0]
  time_naive <-
    as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  #----------------------------------------------------------------------------------
  library(MatrixExtra)
  library(cmfrec)
  X_train <- as.coo.matrix(dat$train.inc)
  # str(X_train)
  X_test <- as.coo.matrix(dat$test.inc)
  # str(X_test)
  
  start_time <- Sys.time()
  model.classic <-
    CMF(
      X_train,
      k = 25,
      lambda = 0.1,
      scale_lam = TRUE,
      verbose = FALSE,
      nthreads = 6
    )
  pred_classic <- predict(model.classic, X_test)
  # print_rmse(X_test, pred_classic, "classic model")
  time_classic = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  start_time <- Sys.time()
  model.baseline <-
    MostPopular(X_train, lambda = 10, scale_lam = FALSE)
  pred_baseline <- predict(model.baseline, X_test)
  # print_rmse(X_test, pred_baseline, "non-personalized model")
  time_baseline <-
    as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  start_time <- Sys.time()
  model.improved <- CMF(
    X_train,
    k = 25,
    lambda = 0.1,
    scale_lam = TRUE,
    add_implicit_features = TRUE,
    w_main = 0.75,
    w_implicit = 0.25,
    use_cg = FALSE,
    niter = 30,
    verbose = FALSE,
    nthreads = 6
  )
  pred_improved <- predict(model.improved, X_test)
  # print_rmse(X_test, pred_improved, "improved classic model")
  time_improved <-
    as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  start_time <- Sys.time()
  model.w.sideinfo <- CMF(
    X_train,
    U = dat$X_r$X,
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
  # print_rmse(X_test, pred_side_info, "model with side info")
  time_sideinfo <-
    as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  #---------------------
  results <- data.frame(
    NonPersonalized = RMSE_error(X_test@x, pred_baseline@x),
    ClassicalModel = RMSE_error(X_test@x, pred_classic@x),
    ClassicPlusImplicit = RMSE_error(X_test@x, pred_improved@x),
    CollectiveModel = RMSE_error(X_test@x, pred_side_info@x),
    CASMAC = RMSE_error(dat$O[dat$W == 0], preds),
    SoftImpute = RMSE_error(dat$O[dat$W == 0], preds_simpute),
    Naive = RMSE_error(dat$O[dat$W == 0], preds_naive)
  )
  results <- as.data.frame(t(results))
  names(results) <- "RMSE"
  results$time <-
    c(
      time_baseline,
      time_classic,
      time_improved,
      time_sideinfo,
      time_casmc,
      time_softImpute,
      time_naive
    )
  
  
  if ("package:cmfrec" %in% search())
    detach("package:cmfrec", unload = TRUE)
  if ("package:MatrixExtra" %in% search())
    detach("package:MatrixExtra", unload = TRUE)
  
  
  results %>%
    kable() %>%
    kable_styling()
}
