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





compare_and_save_with_rep2 <- function(simPar,
                                       num_replications = 500,
                                       n_folds = 5,
                                       new_datfolder = "",
                                       lambda.1_grid = seq(0, 3, length = 20),
                                       lambda.2_grid = seq(.9, 0, length = 20),
                                       alpha_grid = seq(0.992, 1, length = 10),
                                       ncores_mao = 2,
                                       max_cores = 20,
                                       weight_function = Mao_weights$uniform,
                                       n.lambda = 30,
                                       rank.limit = 20,
                                       rank.step = 2,
                                       error_function = error_metric$rmse,
                                       first_seed = NULL,
                                       model_flag = rep(TRUE, 6),
                                       note = "") {
   test_error <<- error_function
   data_dir = paste0("./saved_data/", new_datfolder, "/")
   metrics = c(
      "lambda.beta",
      "lambda.M",
      "time",
      "error.test",
      "corr.test",
      "error.train",
      "error.M",
      "error.beta",
      "rank_M",
      "rank_beta",
      "sparse_in_sparse",
      "nonsparse_in_nonsparse",
      "likelihood_ratio_index",
      "Cox_Snell_R2"
   )
   models = c(
      "SoftImpute",
      "Mao",
      "CASMC-0_Ridge",
      #"CASMC-1_Hard_Rank",
      "CASMC-2_Nuclear",
      "CASMC-3a_Lasso_single_split",
      #"CASMC-3b_Lasso_10-folds",
      "Naive"
   )
   
   model_functions = list(
      SImpute_Sim_Wrapper,
      Mao_Sim_Wrapper,
      CASMC_0_Sim_Wrapper,
      #CASMC_1_Sim_Wrapper,
      CASMC_2_Sim_Wrapper,
      CASMC_3a_Sim_Wrapper,
      #CASMC_3b_Sim_Wrapper,
      Naive_Sim_Wrapper
   )
   # stopifnot(length(model_flag) == length(models))
   # stopifnot(length(model_functions) == length(models))
   # models = models[model_flag]
   # model_functions = model_functions[model_flag]
   
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
            simPar$dim1,
            simPar$dim2,
            r = simPar$nr,
            k = simPar$ncov,
            missing_prob = simPar$missp,
            coll = simPar$coll,
            seed = seed,
            half_discrete = simPar$Disc,
            mar_sparse = simPar$MAR,
            informative_cov_prop = simPar$info,
            mv_beta = TRUE,
            prepare_for_fitting = TRUE
         )
      print(n)
      #----
      # if one model fail, skip to the next, remember to decrease n by 1.
      tryCatch({
         #----------------------------------------------------------------------
         # start fitting the models:
         
         for (i in 1:length(model_functions)) {
            if (i == 1) {
               results <- SImpute_Sim_Wrapper(dat)
               LogLik_SI = results$LogLik
               results = as.numeric(results$results[-1])
            } else
               results = as.numeric(
                  model_functions[[i]](
                     dat = dat,
                     n_folds = n_folds,
                     lambda.1_grid = lambda.1_grid,
                     lambda.2_grid = lambda.2_grid,
                     alpha_grid = alpha_grid,
                     ncores = ncores_mao,
                     max_cores = max_cores,
                     LogLik_SI = LogLik_SI,
                     weight_function = weight_function
                  )[-1]
               )
            new_mean <- update_mean(perf_means[i,], results, n)
            perf_stdev[i,] <-
               update_sse(perf_means[i,], new_mean, perf_stdev[i,],
                          results, n)
            perf_means[i,] <- new_mean
         }
         print(perf_means) |>  round(2)
         # print(perf_stdev)
         
         # I decided to save results after each iteration.
         out_stdev = perf_stdev
         for (i in 1:length(model_functions))
            out_stdev[i,] <- get_sd_from_sse(perf_stdev[i,], n)
         #---------------------
         # we now combine them in a dataframe
         df_means <- as.data.frame(perf_means)
         df_stdev <- as.data.frame(out_stdev)
         # Rename columns
         colnames(df_means) <-
            paste(colnames(df_means), "mean", sep = "_")
         colnames(df_stdev) <-
            paste(colnames(df_stdev), "sd", sep = "_")
         results <- cbind(df_means, df_stdev)
         results$true_rank = dat$rank
         results$dim = paste0("(", simPar$dim1, ",", simPar$dim2, ")")
         results$k = simPar$ncov
         results$B = n
         results$missinginess = simPar$missp
         results$collinearity = simPar$coll
         results$rank_r = simPar$nr
         results$mar_beta = simPar$MAR
         results$half_disc = simPar$Disc
         results$inform_prop = simPar$info
         results$model = models
         
         filename = paste0("Model_Comparison_simulation_replications_",
                           note,
                           ".csv")
         
         write.csv(results,
                   file = paste0(data_dir, filename),
                   row.names = FALSE)
         
      }, error = function(e) {
         print(paste(n, "-", e))
         n = n - 1
      })
      #-------------------------------------
   }
   for (i in 1:length(model_functions))
      perf_stdev[i,] <- get_sd_from_sse(perf_stdev[i,], n)
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
   results$dim = paste0("(", simPar$dim1, ",", simPar$dim2, ")")
   results$k = simPar$ncov
   results$B = n
   results$missinginess = simPar$missp
   results$collinearity = simPar$coll
   results$rank_r = simPar$nr
   results$mar_beta = simPar$MAR
   results$half_disc = simPar$Disc
   results$inform_prop = simPar$info
   results$model = models
   
   #print(results)
   print("Exiting Loop ...")
   
   filename = paste0("Model_Comparison_simulation_replications_",
                     note,
                     ".csv")
   
   write.csv(results,
             file = paste0(data_dir, filename),
             row.names = FALSE)
   print(paste("Results saved to", paste0(data_dir, filename)))
   return(results)
   #----------------------------
}

setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R", local = FALSE)

simPar <- data.frame(
   dim1 = c(193, rep(700,3)),
   dim2 = c(587, rep(800,3)),
   missp = c(.6, rep(0.9,3)),
   coll = FALSE,
   ncov = c(13, 20,20,5),
   nr = 10,
   MAR = c(F, T, F, F),
   Disc = F,
   info = c(.2, .3, .3, 1) 
)

for (i in 2:2) {
   compare_and_save_with_rep2(
      simPar = simPar[i, ],
      num_replications = 30,
      n_folds = 5,
      new_datfolder = "July05",
      lambda.1_grid = seq(1, 0, length = 20),
      lambda.2_grid = seq(.9, 0.1, length = 20),
      alpha_grid = c(1),
      ncores_mao = 1,
      max_cores = 20,
      weight_function = Mao_weights$uniform,
      error_function = error_metric$rmse,
      first_seed =  10,
      model_flag = c(T, T, T, T, T, T, T, T),
      note = paste0("_", i, "_")
   )
    
   
}
