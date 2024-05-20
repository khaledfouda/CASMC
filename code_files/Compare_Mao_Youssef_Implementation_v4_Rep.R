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





compare_and_save_with_rep2 <- function(missingness,
                                      coll = F,
                                      num_replications = 500,
                                      n_folds = 5,
                                      dim = c(400,400),
                                      ncovariates = 10,
                                      lambda.1_grid = seq(0, 3, length = 20),
                                      lambda.2_grid = seq(.9, 0, length = 20),
                                      alpha_grid = seq(0.992, 1, length = 10),
                                      ncores_mao = 2,
                                      max_cores = 20,
                                      weight_function = MaoUniWeights, 
                                      n.lambda = 30,
                                      rank.limit = 20,
                                      rank.step = 2,
                                      error_function = RMSE_error,
                                      first_seed = NULL,
                                      mao_r = ncovariates,
                                      cov_eff = 1,
                                      model_flag = rep(TRUE,6),
                                      note = "") {
   test_error <<- error_function
   data_dir = "./saved_data/"
   stopifnot(missingness %in% c(0, 0.8, 0.9))
   stopifnot(length(dim) == 2)
   stopifnot(length(model_flag) == 6)
   metrics = c(
      "time",
      "lambda.1",
      "lambda.2",
      "error.test",
      "error.all",
      "error.M",
      "error.beta",
      "rank_M",
      "rank_beta"
   )
   models = c(
      "SoftImpute",
      "Mao",
      "CASMC_rank_restriction",
      "CASMC_L2",
      "CASMC_L2_1iter",
      "Naive")
   
   model_functions = list(
      SImpute_Sim_Wrapper,
      Mao_Sim_Wrapper,
      CASMC_rank_Sim_Wrapper,
      function(dat, max_cores,...) CASMC_L2_Sim_Wrapper(dat, max_cores, 300),
      function(dat, max_cores,...) CASMC_L2_Sim_Wrapper(dat, max_cores, 2),
      Naive_Sim_Wrapper
   )
   models = models[model_flag]
   model_functions = model_functions[model_flag]
   
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
         dat <-
         generate_simulation_rows(
            dim[1],
            dim[2],
            ncovariates,
            ncovariates,
            missing_prob = missingness,
            coll = coll,
            seed = seed,
            informative_cov_prop = cov_eff,
            mv_beta = TRUE,
            prepare_for_fitting = TRUE
         )
      print(n)
      #----------------------------------------------------------------------
      # start fitting the models:
      for (i in 1:length(model_functions)) {
         results = as.numeric(
            model_functions[[i]](
               dat = dat,
               n_folds = n_folds,
               lambda.1_grid = lambda.1_grid,
               lambda.2_grid = lambda.2_grid,
               alpha_grid = alpha_grid,
               ncores = ncores_mao,
               max_cores = max_cores,
               weight_function = weight_function
            )[-1]
         )
         new_mean <- update_mean(perf_means[i, ], results, n)
         perf_stdev[i, ] <-
            update_sse(perf_means[i, ], new_mean, perf_stdev[i, ],
                       results, n)
         perf_means[i, ] <- new_mean
      }
      print(perf_means)
      # print(perf_stdev)
      
      # I decided to save results after each iteration.
      out_stdev = perf_stdev
      for (i in 1:length(model_functions))
         out_stdev[i, ] <- get_sd_from_sse(perf_stdev[i, ], n)
      #---------------------
      # we now combine them in a dataframe
      df_means <- as.data.frame(perf_means)
      df_stdev <- as.data.frame(out_stdev)
      # Rename columns
      colnames(df_means) <-
         paste(colnames(df_means), "mean", sep = "_")
      colnames(df_stdev) <- paste(colnames(df_stdev), "sd", sep = "_")
      results <- cbind(df_means, df_stdev)
      results$true_rank = dat$rank
      results$dim = paste0("(",dim[1],",",dim[2],")") 
      results$k = ncovariates
      results$B = n
      results$missinginess = missingness
      results$collinearity = coll
      results$model = models
      
      filename = paste0(
         "Compare_MC_Models_Youssef_Simulation_with_replications_",
         note,
         round(missingness * 100),
         "_coll_",
         coll,
         "_dim_",
         dim[1],
         "x",
         dim[2],
         ".csv"
      )
      
      write.csv(results,
                file = paste0(data_dir, filename),
                row.names = FALSE)
      
      
      #-------------------------------------
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
   results$true_rank = dat$rank
   results$dim = paste0("(",dim[1],",",dim[2],")")
   results$k = ncovariates
   results$B = num_replications
   results$missinginess = missingness
   results$collinearity = coll
   results$model = models
   
   #print(results)
   print("Exiting Loop ...")
   
      filename = paste0(
         "Compare_MC_Models_Youssef_Simulation_with_replications_",
         note,
         round(missingness * 100),
         "_coll_",
         coll,
         "_dim_",
         dim[1],
         "x",
         dim[2],
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


for (hparams in list(#c(0.9, 400,500),
                     c(0.9, 700, 800)
                     # c(0.9, TRUE),
                     # c(0.8, TRUE)
                     #c(0.0, FALSE)
                     )) {
   compare_and_save_with_rep2(
      hparams[1],
      FALSE,
      num_replications = 50,
      dim = c(hparams[2],hparams[3]),
      lambda.1_grid = seq(1, 0, length = 20),
      lambda.2_grid = seq(.9, 0.1, length = 20),
      alpha_grid = c(1),
      ncores_mao = 1,
      max_cores = 20,
      ncovariates = 10,
      mao_r = 10,
      weight_function = Mao_weights$uniform,
      error_function = error_metric$rmse, 
      cov_eff = 0.7,
      first_seed =  10,
      model_flag = c(T,T,T,T,T,T),
      note = "_new47_"
   )
}
