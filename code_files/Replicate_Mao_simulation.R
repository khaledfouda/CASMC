compare_Mao_mod <-
  function(gen.dat,
           lambda.1_grid,
           lambda.2_grid,
           alpha_grid,
           numCores = 1,
           n_folds = 5) {
    start_time = Sys.time()
    cv.out <- Mao.cv(
      gen.dat$A,
      gen.dat$X,
      gen.dat$Y,
      gen.dat$W,
      n_folds = n_folds,
      lambda.1_grid = lambda.1_grid,
      lambda.2_grid = lambda.2_grid,
      alpha_grid = alpha_grid,
      numCores = ncores,
      n1n2_optimized = FALSE,
      theta_estimator = theta_default
    )
    mao.out <-
      Mao.fit(
        gen.dat$Y,
        gen.dat$X,
        gen.dat$W,
        cv.out$best_parameters$lambda.1,
        cv.out$best_parameters$lambda.2,
        cv.out$best_parameters$alpha,
        theta_estimator = theta_default,
        n1n2_optimized = FALSE
      )
    results = list()
    results$error.beta = RMSE_error(mao.out$beta_hat, gen.dat$beta)
    results$error.B = RMSE_error(mao.out$B_hat, gen.dat$B)
    results$error.all = RMSE_error(mao.out$A_hat, gen.dat$A)
    results$error.test = mao_error(mao.out$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
                                                                              0])
    results$rank = mao.out$rank
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results
  }


compare_softImpute_orig_mod <- function(gen.dat, valid.dat) {
  start_time = Sys.time()
  sout <- simpute.orig(
    valid.dat$Y_train,
    valid.dat$W_valid,
    gen.dat$Y,
    trace = FALSE,
    rank.limit = 30,
    print.best = FALSE,
    rank.step = 4
  )
  results = list()
  results$error.beta = NA
  results$error.B = NA
  results$error.all = RMSE_error(sout$A_hat, gen.dat$A)
  results$error.test = mao_error(sout$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
                                                                         0])
  results$rank = sout$rank_A
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results
}



compare_softImpute_cov_mod <- function(gen.dat, valid.dat) {
  start_time = Sys.time()
  sout <-
    simpute.cov.cv(
      valid.dat$Y_train,
      gen.dat$X,
      valid.dat$W_valid,
      valid.dat$Y_valid,
      trace = F,
      rank.limit = 30,
      quiet = TRUE,
      print.best = FALSE,
      rank.step = 4,
      type = "als",
      lambda1 = 0,
      tol = 2
    )
  results = list()
  results$error.beta = RMSE_error(sout$beta_hat, gen.dat$beta)
  results$error.B = RMSE_error(sout$B_hat, gen.dat$B)
  results$error.all = RMSE_error(sout$A_hat, gen.dat$A)
  results$error.test = mao_error(sout$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
                                                                         0])
  results$rank = sout$rank_A
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results
}


compare_softImpute_splr_mod <-
  function(gen.dat,
           valid.dat,
           splr.dat,
           rank.step,
           rank.limit,
           n.lambda) {
    start_time = Sys.time()
    
    
    best_fit = simpute.cov.cv_splr(
      valid.dat$Y_train,
      splr.dat,
      valid.dat$Y_valid,
      valid.dat$W_valid,
      gen.dat$Y,
      trace = F,
      thresh = 1e-6,
      n.lambda = n.lambda,
      rank.limit = rank.limit,
      maxit = 200,
      rank.step = rank.step,
      print.best = FALSE
    )
    fit1 = best_fit$fit1
    fit2 = best_fit$fit2
    sout = best_fit
    # get estimates and validate
    sout$B_hat = fit1$u %*% (fit1$d * t(fit1$v))
    sout$beta_hat =  fit2$u %*% (fit2$d * t(fit2$v))
    sout$A_hat = sout$B_hat + splr.dat$X %*% t(sout$beta_hat)
    
    results = list()
    results$error.beta = RMSE_error(t(sout$beta_hat), gen.dat$beta)
    results$error.B = RMSE_error(sout$B_hat, gen.dat$B)
    results$error.all = RMSE_error(sout$A_hat, gen.dat$A)
    results$error.test = mao_error(sout$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
                                                                           0])
    results$rank = qr(sout$A_hat)$rank
    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    results
  }





replicate_mao_sim <- function(n_folds = 3,
                              n_reps = 500,
                              dim = 400,
                              ncovariates = 20,
                              mao_r = 10,
                              lambda.1_grid = seq(0, 2, length = 20),
                              lambda.2_grid = seq(.9, 0, length = 20),
                              alpha_grid = seq(0.992, 1, length = 10),
                              ncores = 1,
                              n.lambda = 30,
                              rank.limit = 20,
                              rank.step = 2,
                              error_function = RMSE_error,
                              save_to_disk = FALSE,
                              seed = NULL,
                              print_every = 10) {
  print(paste("Running on", ncores, "core(s)."))
  test_error <<- error_function
  data_dir = "./saved_data/"
  results1 = results2 = results3 = results4 = data.frame(
    error_beta = rep(NA, n_reps),
    error_B = NA,
    error_A = NA,
    error_test = NA,
    rank = NA,
    time = NA
  )
  
  seed = ifelse(is.null(seed), 0, seed)
  for (i in 1:n_reps) {
    seed = seed + i
    #final.results <- foreach(i = 1:length(dim), .combine='rbind') %do%  {
    set.seed(seed)
    gen.dat <-
      generate_simulation_data_mao(
        n1 = dim,
        n2 = dim,
        m = ncovariates,
        r = mao_r,
        seed = seed
      )
    #-------------------------------------------------------------------------------------
    # validation set to be used for the next two models
    valid.dat = list()
    valid.dat$W_valid <-
      matrix.split.train.test(gen.dat$W, testp = 0.2, seed = seed)
    valid.dat$Y_train <- gen.dat$Y * valid.dat$W_valid
    valid.dat$Y_valid <- gen.dat$Y[valid.dat$W_valid == 0]
    #---------------------------------------------------------
    # SPLR data
    splr.dat = reduced_hat_decomp(gen.dat$X, 1e-2)
    gen.dat$X <- splr.dat$X
    if (i %% print_every == 0)
      print(i)
    #----------------------------------------------------------------------
    # fit 1. Mao
    
    results1[i,] <-  compare_Mao_mod(gen.dat,
                                     lambda.1_grid,
                                     lambda.2_grid,
                                     alpha_grid,
                                     1,
                                     n_folds)
    
    
    #----------------------------------------------------------
    # soft Impute model without covariates
    results2[i,] <- compare_softImpute_orig_mod(gen.dat, valid.dat)
    #----------------------------------------------------------------------------
    # soft Impute model with covariates
    results3[i,] <- compare_softImpute_cov_mod(gen.dat, valid.dat)
    #--------------------------------------------------------------------------------
    results4[i,] <- compare_softImpute_splr_mod(gen.dat,
                                                valid.dat,
                                                splr.dat,
                                                rank.step,
                                                rank.limit,
                                                n.lambda)
    #---------------------------------------------------------------------------
  }
  print("Exiting Loop ...")
  
  
  results = rbind(
    summarise_all(results1, list(
      mean = ~ mean(., na.rm = TRUE),
      sd = ~ sd(., na.rm = TRUE)
    )),
    summarise_all(results2, list(
      mean = ~ mean(., na.rm = TRUE),
      sd = ~ sd(., na.rm = TRUE)
    )),
    summarise_all(results3, list(
      mean = ~ mean(., na.rm = TRUE),
      sd = ~ sd(., na.rm = TRUE)
    )),
    summarise_all(results4, list(
      mean = ~ mean(., na.rm = TRUE),
      sd = ~ sd(., na.rm = TRUE)
    ))
  )
  
  results$dim = dim
  results$ncovariates = ncovariates
  results$n_folds = n_folds
  results$n_reps = n_reps
  results$true_rank = gen.dat$rank
  
  if (save_to_disk) {
    filename = paste0("Replicate_Mao_Simulation_",
                      n_reps,
                      "_",
                      n_folds,
                      "_",
                      dim,
                      ".csv")
    write.csv(results,
              file = paste0(data_dir, filename),
              row.names = FALSE)
    print(paste("Results saved to", paste0(data_dir, filename)))
    
  }
  print(results)
  return(results)
  #----------------------------
}
#--------------------------
#
#
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R", local = FALSE)


for(dim in c(400,600,800))
results = replicate_mao_sim(
  n_folds = 3,
  n_reps = 200,
  dim = dim,
  save_to_disk = TRUE,
  print_every = 10,
  error_function = RMSE_error 
)
