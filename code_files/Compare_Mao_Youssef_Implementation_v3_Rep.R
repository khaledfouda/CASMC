require(foreach)
require(doParallel)

update_mean <- function(current_mean, new_data, n) {
   current_mean + (new_data - current_mean) / n
}
update_sse <-
   function(current_mean,
            new_mean,
            current_sse,
            new_data,
            n) {
      current_sse + (new_data - current_mean) * (new_data - new_mean)
   }
get_sd_from_sse <- function(current_sse, n) {
   sqrt(current_sse / (n - 1))
}





compare_and_save_with_rep <- function(missingness,
                                      coll = TRUE,
                                      num_replications = 500,
                                      n_folds = 5,
                                      dim = 400,
                                      ncovariates = 10,
                                      lambda.1_grid = seq(0, 3, length = 20),
                                      lambda.2_grid = seq(.9, 0, length = 20),
                                      alpha_grid = seq(0.992, 1, length = 10),
                                      ncores_mao = 2,
                                      weight_function = MaoUniWeights,
                                      n.lambda = 30,
                                      rank.limit = 20,
                                      rank.step = 2,
                                      error_function = RMSE_error,
                                      first_seed = NULL,
                                      mao_r = ncovariates,
                                      cov_eff = TRUE,
                                      note = "") {
   test_error <<- error_function
   data_dir = "./saved_data/"
   stopifnot(missingness %in% c(0, 0.8, 0.9))
   stopifnot(length(dim) == 1)
   metrics = c(
      "time",
      "alpha",
      "lambda.1",
      "lambda.2",
      "error.test",
      "error.all",
      "error.M",
      "error.beta",
      "rank"
   )
   models = c(
      "Mao",
      "SoftImpute",
      "CAMC_holdout",
      "CAMC_kfold",
      "CASMC_holdout",
      "CASMC_kfold",
      "Naive")
   
   model_functions = list(
      compare_Mao,
      compare_softImpute,
      compare_CAMC_holdout,
      compare_CAMC_kfold,
      compare_CASMC_holdout,
      compare_CASMC_kfold,
      compare_naive
   )
   perf_means <- perf_stdev <-
      matrix(0,
             length(models),
             length(metrics),
             dimnames = list(models, metrics))
   seed = ifelse(is.null(first_seed), 0, first_seed)
   
   for (n in 1:num_replications) {
      # set a different seed at each step
      seed = seed + n
      set.seed(seed)
      #----------------------------------
      # results <- data.frame()
      
      #---- initialize the data
      if (missingness == 0) {
         gen.dat <-
            generate_simulation_data_mao(
               n1 = dim,
               n2 = dim,
               m = ncovariates,
               r = mao_r,
               seed = seed,
               cov_eff = cov_eff
            )
      } else
         gen.dat <-
         generate_simulation_data_ysf(
            2,
            dim,
            dim,
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
      #gen.dat$X <- splr.dat$X
      print(n)
      #----------------------------------------------------------------------
      # start fitting the models:
      for (i in 1:length(model_functions)) {
         results = as.numeric(
            model_functions[[i]](
               gen.dat = gen.dat,
               valid.dat = valid.dat,
               splr.dat = splr.dat,
               n_folds = n_folds,
               rank.step = rank.step,
               rank.limit = rank.limit,
               n.lambda = n.lambda,
               lambda.1_grid = lambda.1_grid,
               lambda.2_grid = lambda.2_grid,
               alpha_grid = alpha_grid,
               ncores = ncores_mao,
               weight_function = weight_function
            )[-1]
         )
         new_mean <- update_mean(perf_means[i, ], results, n)
         perf_stdev[i, ] <-
            update_sse(perf_means[i, ], new_mean, perf_stdev[i, ],
                       results, n)
         perf_means[i, ] <- new_mean
      }
      # print(perf_means)
      # print(perf_stdev)
   }
   for (i in 1:length(model_functions))
      perf_stdev[i, ] <- get_sd_from_sse(perf_stdev[i, ], n)
   #---------------------
   # we now combine them in a dataframe
   df_means <- as.data.frame(perf_means)
   df_stdev <- as.data.frame(perf_stdev)
   # Rename columns
   colnames(df_means) <-
      paste(colnames(df_means), "mean", sep = "_")
   colnames(df_stdev) <- paste(colnames(df_stdev), "sd", sep = "_")
   results <- cbind(df_means, df_stdev)
   results$true_rank = gen.dat$rank
   results$dim = dim
   results$k = ncovariates
   results$B = num_replications
   results$missinginess = missingness
   results$collinearity = coll
   
   print(results)
   print("Exiting Loop ...")
   
   
   if (missingness == 0) {
      filename = paste0("Compare_MC_Models_Mao_Simulation_with_replications",
                        note,
                        ".csv")
   } else
      filename = paste0(
         "Compare_MC_Models_Youssef_Simulation_with_replications_",
         note,
         round(missingness * 100),
         "_coll_",
         coll,
         "_dim_",
         dim,
         ".csv"
      )
   
   write.csv(results,
             file = paste0(data_dir, filename),
             row.names = FALSE)
   print(paste("Results saved to", paste0(data_dir, filename)))
   return(results)
   #----------------------------
}

setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R", local = FALSE)


for (hparams in list(c(0.9, TRUE),
                     c(0.9, FALSE),
                     c(0, FALSE),
                     c(0.8, TRUE),
                     c(0.8, FALSE))) {
   compare_and_save_with_rep(
      hparams[1],
      hparams[2],
      num_replications = 50,
      dim = 600,
      lambda.1_grid = seq(2, 0, length = 20),
      lambda.2_grid = seq(.9, 0, length = 20),
      alpha_grid = seq(0.992, 1, length = 10),
      ncores_mao = 1,
      ncovariates = 10,
      mao_r = 10,
      error_function = RMSE_error,
      cov_eff = TRUE,
      note = ""
   )
}
