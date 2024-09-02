

setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
#---------------
print_rmse <- function(X_test, X_hat, model_name) {
  rmse <- sqrt(mean((X_test@x - X_hat@x) ^ 2))
  cat(sprintf("RMSE for %s is: %.4f\n", model_name, rmse))
}

RMSE_error2 <- function(preds, orig, glob_mean, glob_sd) {
  error_metric$rmse((preds * glob_sd) + glob_mean, (orig * glob_sd) + glob_mean)
}

subs <- c("a", "b")
i = 1
all_res <- list()
for (i in 1:2) {
  source("./code_files/import_lib.R")
  source("./MovieLens/load_data_100k.R")
  dat <-
    load_movielens_100k(subs[i], remove_bad_movies = T, scale = F)
  dat$test.ind = cbind(dat$test.df$user_id, dat$test.df$movie_id)
  
  if (i == 2) {
    best_fit = CASMC_cv_holdout_with_r(
      dat$valid$train.inc,
      dat$X_r,
      dat$valid$test,
      dat$valid$W,
      dat$train.inc,
      #r = dat$X_r$rank,
      r_min = 0,
      r_max =  15,
      track_r = T,
      trace = T,
      thresh = 1e-6,
      n.lambda = 30,
      lambda.factor = 0.9,
      rank.init = 1,
      rank.step = 2,
      rank.limit = 30,
      pct = 0.90,
      seed = 2023,
      error_function = error_metric$rmse
    )
  } else{
    best_fit = CASMC_cv_holdout_with_r(
      dat$valid$train.inc,
      dat$X_r,
      dat$valid$test,
      dat$valid$W,
      dat$train.inc,
      # r = dat$X_r$rank,
      r_min = 0,
      r_max =  15,
      track_r = T,
      trace = T,
      thresh = 1e-6,
      n.lambda = 10,
      rank.init = 5,
      rank.step = 1,
      rank.limit = 30,
      pct = 0.90,
      error_function = error_metric$rmse
    )
  }
  
  
  fit1 = best_fit$fit
  beta =  fit1$Beta$u %*% (fit1$Beta$d * t(fit1$Beta$v))
  M = fit1$u %*% (fit1$d * t(fit1$v))
  A = M + dat$X_r$X %*% t(beta)
  #A = revertBiScaledMatrix(as.matrix(A), dat$biScale)
  #qr(A)$rank
  preds <- A[dat$test.ind]
  error_metric$rmse((preds * dat$glob_sd) + dat$glob_mean,
             (dat$test.df$rating * dat$glob_sd) + dat$glob_mean
  ) %>% print()
  error_metric$rmse(preds, dat$test.df$rating) %>% print()
  
  
  best_fit$r %>% print()
  #------------------------------------------------------------
  # soft-impute
  start_time = Sys.time()
  sout <- simpute.cv(
    dat$valid$train,
    dat$train.mat,
    trace = FALSE,
    rank.limit = 30,
    print.best = FALSE,
    rank.step = 2,
    test_error = error_metric$rmse
  )
  #sout$estimates <-  revertBiScaledMatrix(as.matrix(sout$estimates), dat$biScale)
  preds_simpute <- sout$estimates[dat$test.ind]
  
  time_softImpute = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  #---------------------------------------------------------------------------------
  # Naive
  start_time = Sys.time()
  estimates = as.matrix(naive_MC(dat$train.mat))
  #estimates =  revertBiScaledMatrix(as.matrix(estimates), dat$biScale)
  preds_naive = estimates[dat$test.ind]
  time_naive <-
    as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  #----------------------------------------------------------------------------------
  
  
  
  #---------------------
  library(cmfrec)
  library(Matrix)
  library(MatrixExtra)
  
  # dat <-
  #  load_movielens_100k("a", remove_bad_movies = TRUE, scale = FALSE)
  
  
  
  
  X_train <- as.coo.matrix(dat$train.inc)
  str(X_train)
  X_test <- as.coo.matrix(dat$test.inc)
  str(X_test)
  
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
  print_rmse(X_test, pred_classic, "classic model")
  
  model.baseline <-
    MostPopular(X_train, lambda = 10, scale_lam = FALSE)
  pred_baseline <- predict(model.baseline, X_test)
  print_rmse(X_test, pred_baseline, "non-personalized model")
  
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
  print_rmse(X_test, pred_improved, "improved classic model")
  
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
  print_rmse(X_test, pred_side_info, "model with side info")
  detach("package:cmfrec", unload = TRUE)
  detach("package:MatrixExtra", unload = TRUE)
  
  #------------------------------------------------
  results <- data.frame(
    NonPersonalized = RMSE_error2(X_test@x, pred_baseline@x, dat$glob_mean, dat$glob_sd),
    ClassicalModel = RMSE_error2(X_test@x, pred_classic@x, dat$glob_mean, dat$glob_sd),
    ClassicPlusImplicit = RMSE_error2(X_test@x, pred_improved@x, dat$glob_mean, dat$glob_sd),
    CollectiveModel = RMSE_error2(X_test@x, pred_side_info@x, dat$glob_mean, dat$glob_sd),
    CASMAC = RMSE_error2(dat$test.df$rating, preds, dat$glob_mean, dat$glob_sd),
    SoftImpute = RMSE_error2(
      dat$test.df$rating,
      preds_simpute,
      dat$glob_mean,
      dat$glob_sd
    ),
    Naive = RMSE_error2(dat$test.df$rating, preds_naive, dat$glob_mean, dat$glob_sd)
  )
  results <- as.data.frame(t(results))
  names(results) <- "RMSE"
  
  rmse_naive = RMSE_error2(dat$test.df$rating, preds_naive, dat$glob_mean, dat$glob_sd)
  
  improvements <-
    sapply(results, function(x)
      (rmse_naive - x) / rmse_naive * 100) %>% round(1)
  results$improvement = (improvements)
  all_res[[i]] <- results
}


all_res[[1]] %>%
  arrange(desc(improvement)) %>% 
  kable() %>%
  kable_styling()



all_res[[2]] %>%
  arrange(desc(improvement)) %>% 
  kable() %>%
  kable_styling()
