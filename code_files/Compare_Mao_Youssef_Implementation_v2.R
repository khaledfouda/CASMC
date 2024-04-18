require(foreach)
require(doParallel)



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
   stopifnot(length(model_mask) == 7)
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
               n_folds,
               weight_function = MaoUniWeights
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

alpha_grid = seq(0.992, 1, length = 5)
lambda.1_grid = seq(20, 50, length = 20)
lambda.2_grid = seq(.9, 0, length = 20)
ncores = 1
error_function <- RMSE_error
model_mask <- rep(T, 7)
model_mask[c(4)] <- F
mao_r <- 10
ncovariates <- 5
cov_eff = F
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
