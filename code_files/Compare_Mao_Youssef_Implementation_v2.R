require(foreach)
require(doParallel)

compare_Mao <-
   function(gen.dat,
            lambda.1_grid,
            lambda.2_grid,
            alpha_grid,
            ncores = 1,
            n_folds = 5,
            weight_function = MaoBinomalWeights) {
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
            n1n2_optimized = FALSE
         )
      
      results = list(model = "Mao")
      results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
      results$alpha = cv.out$best_parameters$alpha
      results$lambda.1 = cv.out$best_parameters$lambda.1
      results$lambda.2 = cv.out$best_parameters$lambda.2
      results$error.test = test_error(mao.out$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
                                                                                  0])
      results$error.all = test_error(mao.out$A_hat, gen.dat$A)
      results$error.B = test_error(mao.out$B_hat, gen.dat$B)
      results$error.beta = test_error(mao.out$beta_hat, gen.dat$beta)
      results$rank = mao.out$rank
      results
   }

compare_softImpute <- function(gen.dat, valid.dat) {
   start_time = Sys.time()
   sout <- simpute.cv(
      valid.dat$Y_train,
      valid.dat$W_valid,
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
   results$error.test = test_error(sout$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
                                                                            0])
   results$error.all = test_error(sout$A_hat, gen.dat$A)
   results$error.B = NA
   results$error.beta = NA
   results$rank = sout$rank_A
   results
}

# compare_softImpute_cov <- function(gen.dat, valid.dat) {
#    start_time = Sys.time()
#    sout <-
#       simpute.cov.cv(
#          valid.dat$Y_train,
#          gen.dat$X,
#          valid.dat$W_valid,
#          valid.dat$Y_valid,
#          trace = F,
#          rank.limit = 30,
#          quiet = TRUE,
#          print.best = FALSE,
#          rank.step = 4,
#          type = "als",
#          lambda1 = 0,
#          tol = 2
#       )
#    results = list(model = "SoftImpute_Cov")
#    results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
#    results$alpha = NA
#    results$lambda.1 = 0
#    results$lambda.2 = sout$lambda
#    results$error.test = test_error(sout$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
#                                                                             0])
#    results$error.all = test_error(sout$A_hat, gen.dat$A)
#    results$error.B = test_error(sout$B_hat, gen.dat$B)
#    results$error.beta = test_error(sout$beta_hat, gen.dat$beta)
#    results$rank = sout$rank_A
#    results
# }

compare_CAMC_holdout <-
   function(gen.dat,
            valid.dat,
            lambda.1_grid,
            rank.step,
            rank.limit,
            n.lambda) {
      start_time = Sys.time()
      
      fiti <- CAMC_cv_holdout(valid.dat$Y_train, gen.dat$X, valid.dat$W_valid,
                              valid.dat$Y_valid, trace = FALSE, rank.limit = rank.limit,
                              print.best = FALSE, rank.step  = rank.step, n.lambda = n.lambda,
                              type = "als", quiet = TRUE, tol = 2, lambda.1_grid = lambda.1_grid)
      sout <- fiti$fit2
         
      
      results = list(model = "CAMC_holdout")
      results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
      results$alpha = NA
      results$lambda.1 = sout$lambda1
      results$lambda.2 = sout$lambda2
      results$error.test = test_error(sout$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
                                                                               0])
      results$error.all = test_error(sout$A_hat, gen.dat$A)
      results$error.beta = test_error(sout$beta_hat, gen.dat$beta)
      results$error.B = test_error(sout$B_hat, gen.dat$B)
      results$rank = sout$rank_A
      results
   }

compare_CAMC_kfold <-
   function(gen.dat,
            lambda.1_grid,
            n_folds,
            rank.step,
            rank.limit,
            n.lambda) {
      start_time = Sys.time()
      fiti <- CAMC_cv_kfold(gen.dat$Y,
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
                            tol = 2)
      sout <- fiti$fit2
      results = list(model = "CAMC_kfold")
      results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
      results$alpha = NA
      results$lambda.1 = sout$lambda1
      results$lambda.2 = sout$lambda2
      results$error.test = test_error(sout$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
                                                                               0])
      results$error.all = test_error(sout$A_hat, gen.dat$A)
      results$error.beta = test_error(sout$beta_hat, gen.dat$beta)
      results$error.B = test_error(sout$B_hat, gen.dat$B)
      results$rank = sout$rank_A
      results
   }




compare_CASMC_holdout <-
   function(gen.dat,
            valid.dat,
            splr.dat,
            rank.step,
            rank.limit,
            n.lambda) {
      start_time = Sys.time()
      
      
      best_fit = CASMC_cv_holdout_v2(
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
      
      
      fit1 = best_fit$fit
      #fit2 = best_fit$fit2
      sout = best_fit
      # get estimates and validate
      sout$B_hat = fit1$u %*% (fit1$d * t(fit1$v))
      sout$beta_hat =  fit1$beta#fit2$u %*% (fit2$d * t(fit2$v))
      sout$A_hat = sout$B_hat + splr.dat$X %*% t(sout$beta_hat)
      
      results = list(model = "CASMC_holdout")
      results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
      results$alpha = NA
      results$lambda.1 = NA
      results$lambda.2 = sout$lambda
      results$error.test = test_error(sout$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
                                                                               0])
      results$error.all = test_error(sout$A_hat, gen.dat$A)
      results$error.beta = test_error(t(sout$beta_hat), gen.dat$beta)
      results$error.B = test_error(sout$B_hat, gen.dat$B)
      results$rank = qr(sout$A_hat)$rank
      results
   }

compare_CASMC_kfold <-
   function(gen.dat,
            splr.dat,
            n_folds,
            rank.step,
            rank.limit,
            n.lambda) {
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
      #fit2 = best_fit$fit2
      sout = best_fit
      # get estimates and validate
      sout$B_hat = fit1$u %*% (fit1$d * t(fit1$v))
      sout$beta_hat =  fit1$beta#fit2$u %*% (fit2$d * t(fit2$v))
      sout$A_hat = sout$B_hat + splr.dat$X %*% t(sout$beta_hat)
      
      
      results = list(model = "CASMC_kfold")
      results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
      results$alpha = NA
      results$lambda.1 = NA
      results$lambda.2 = sout$lambda
      results$error.test = test_error(sout$A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W ==
                                                                               0])
      results$error.all = test_error(sout$A_hat, gen.dat$A)
      results$error.beta = test_error(t(sout$beta_hat), gen.dat$beta)
      results$error.B = test_error(sout$B_hat, gen.dat$B)
      results$rank = qr(sout$A_hat)$rank
      results
   }

compare_naive <- function(gen.dat) {
   start_time = Sys.time()
   A_hat = naive_MC(gen.dat$Y)
   results = list(model = "Naive")
   results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
   results$alpha = NA
   results$lambda.1 = NA
   results$lambda.2 = NA
   results$error.test = test_error(A_hat[gen.dat$W == 0], gen.dat$A[gen.dat$W == 0])
   results$error.all = test_error(A_hat, gen.dat$A)
   results$error.B = NA
   results$error.beta = NA
   results$rank = qr(A_hat)$rank
   results
}


compare_and_save <- function(missingness,
                             coll = TRUE,
                             n_folds = 5,
                             dim = seq(400, 1000, 200),
                             ncovariates = 10,
                             lambda.1_grid = seq(0, 3, length = 20),
                             lambda.2_grid = seq(.9, 0, length = 20),
                             alpha_grid = seq(0.992, 1, length = 10),
                             ncores = 1,
                             n.lambda = 30,
                             rank.limit = 20,
                             rank.step = 2,
                             error_function = RMSE_error,
                             seed = NULL,
                             model_mask = rep(TRUE, 8),
                             mao_r = ncovariates,
                             cov_eff = TRUE,
                             note = "") {
   ncores = min(ncores, length(dim))
   stopifnot(length(model_mask)==7)
   if (ncores > 1) {
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
   } else
      registerDoSEQ()
   
   print(paste("Running on", ncores, "core(s)."))
   test_error <<- error_function
   data_dir = "./saved_data/"
   stopifnot(missingness %in% c(0, 0.8, 0.9))
   final.results = data.frame()
   
   for (i in 1:length(dim)) {
      #final.results <- foreach(i = 1:length(dim), .combine='rbind') %do%  {
      if (!is.null(seed))
         set.seed(seed)
      results <- data.frame()
      if (missingness == 0) {
         gen.dat <-
            generate_simulation_data_mao(
               n1 = dim[i],
               n2 = dim[i],
               m = ncovariates,
               r = mao_r,
               seed = seed,
               cov_eff = cov_eff
            )
      } else
         gen.dat <-
         generate_simulation_data_ysf(
            2,
            dim[i],
            dim[i],
            ncovariates,
            ncovariates,
            missing_prob = missingness,
            coll = coll,
            seed = seed,
            cov_eff = cov_eff
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
      print(i)
      #----------------------------------------------------------------------
      # fit 1. Mao
      cat(" - M1 - ")
      if (model_mask[1])
         results = rbind(
            results,
            compare_Mao(
               gen.dat,
               lambda.1_grid,
               lambda.2_grid,
               alpha_grid,
               1,
               n_folds
            )
         )
      #----------------------------------------------------------
      # soft Impute model without covariates
      cat("M2 - ")
      if (model_mask[2])
         results = rbind(results, compare_softImpute(gen.dat, valid.dat))
      #----------------------------------------------------------------------------
      # Soft Impute with Covariates and With L2 regularization on the covariates
      cat("M3 - ")
      if (model_mask[3])
         results = rbind(
            results,
            compare_CAMC_holdout(
               gen.dat,
               valid.dat,
               lambda.1_grid,
               rank.step,
               rank.limit,
               n.lambda
            )
         )
      #-------------------------------------------------------------------------------------
      # Soft Impute with Covariates and With L2 regularization on the covariates and K-fold cross-validation
      cat("M4 - ")
      if (model_mask[4])
         results = rbind(
            results,
            compare_CAMC_kfold(
               gen.dat,
               lambda.1_grid,
               n_folds,
               rank.step,
               rank.limit,
               n.lambda
            )
         )
      #--------------------------------------------------------------------------------
      cat("M5 - ")
      if (model_mask[5])
         results = rbind(
            results,
            compare_CASMC_holdout(
               gen.dat,
               valid.dat,
               splr.dat,
               rank.step,
               rank.limit,
               n.lambda
            )
         )
      #--------------------------------------------------------------------------------
      cat("M6 - ")
      if (model_mask[6])
         results = rbind(
            results,
            compare_CASMC_kfold(
               gen.dat,
               splr.dat,
               n_folds,
               rank.step,
               rank.limit,
               n.lambda
            )
         )
      #--------------------------------------------------------------------------------
      cat("M7 - ")
      if (model_mask[7])
         results = rbind(results, compare_naive(gen.dat))
      cat("Done.\n")
      #--------------------------------------------------------------------------------
      results$true_rank = gen.dat$rank
      results$dim = dim[i]
      results$k = ncovariates
      #---------------------------------------------------------------------------
      # saving plots to disk
      print(results)
      #return(results)
      final.results = rbind(final.results, results)
   }
   print("Exiting Loop ...")
   if (ncores > 1)
      stopCluster(cl)
   #stopImplicitCluster()
   
   final.results$missinginess = missingness
   final.results$collinearity = coll
   
   if (missingness == 0) {
      filename = paste0("Compare_MC_Models_Mao_Simulation", note, ".csv")
   } else
      filename = paste0(
         "Compare_MC_Models_Youssef_Simulation_",
         note,
         round(missingness * 100),
         "_coll_",
         coll,
         ".csv"
      )
   
   write.csv(final.results,
             file = paste0(data_dir, filename),
             row.names = FALSE)
   print(paste("Results saved to", paste0(data_dir, filename)))
   return(final.results)
   #----------------------------
} 

setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R", local = FALSE)

alpha_grid = seq(0.992, 1, length = 10)
lambda.1_grid = seq(0, 2, length = 20)
lambda.2_grid = seq(.9, 0, length = 20)
ncores = 1
error_function <- RMSE_error
model_mask <- rep(T, 7)
model_mask[c(1, 2,6)] <- TRUE
mao_r <- 10
ncovariates <- 10
cov_eff = FALSE
note = "_no_cov_"


compare_and_save(
   0.8,
   FALSE,
   lambda.1_grid = lambda.1_grid,
   lambda.2_grid = lambda.2_grid,
   alpha_grid = alpha_grid,
   ncores = ncores,
   model_mask = model_mask,
   ncovariates = ncovariates,
   mao_r = mao_r,
   error_function = error_function,
   cov_eff = cov_eff,
   note = note
)
compare_and_save(
   0.8,
   TRUE,
   lambda.1_grid = lambda.1_grid,
   lambda.2_grid = lambda.2_grid,
   alpha_grid = alpha_grid,
   ncores = ncores,
   model_mask = model_mask,
   ncovariates = ncovariates,
   mao_r = mao_r,
   error_function = error_function,
   cov_eff = cov_eff,
   note = note
)
compare_and_save(
   0.9,
   FALSE,
   lambda.1_grid = lambda.1_grid,
   lambda.2_grid = lambda.2_grid,
   alpha_grid = alpha_grid,
   ncores = ncores,
   model_mask = model_mask,
   ncovariates = ncovariates,
   mao_r = mao_r,
   error_function = error_function,
   cov_eff = cov_eff,
   note = note
)
compare_and_save(
   0.9,
   TRUE,
   lambda.1_grid = lambda.1_grid,
   lambda.2_grid = lambda.2_grid,
   alpha_grid = alpha_grid,
   ncores = ncores,
   model_mask = model_mask,
   ncovariates = ncovariates,
   mao_r = mao_r,
   error_function = error_function,
   cov_eff = cov_eff,
   note = note
)
compare_and_save(
   0,
   FALSE,
   lambda.1_grid = lambda.1_grid,
   lambda.2_grid = lambda.2_grid,
   alpha_grid = alpha_grid,
   ncores = ncores,
   model_mask = model_mask,
   ncovariates = ncovariates,
   mao_r = mao_r,
   error_function = error_function,
   cov_eff = cov_eff,
   note = note
)
#------------------------------------------------------------------------------------
