



setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
source("./code_files/import_lib.R")
source("./BIXI/data-raw/bixi_data.R")
source("./BIXI/model_comparison/fit_wrappers_bixi.R")

num_replications = 5
spatial = FALSE
if(spatial){
dat <-
 load_bixi_dat(transpose = F, scale_response = T, scale_covariates = F,
               testp = 0.2, validp = 0.2, seed=2023)$model
dat$X <- dat$spatial_simple
dat$X <- dat$X |> #[, c(1, 2, 5)] |>
  scalers("minmax")  
}else{
dat <-
 load_bixi_dat(
  transpose = T,
  scale_response = T,
  scale_covariates = F,
  testp = 0.2,
  validp = 0.2,
  seed = 2023
 )$model

dat$X <- dat$X[, c(1, 2, 5)] |>
  scalers("minmax")
}

# 


dat$masks$tr_val = (dat$masks$obs == 1) & (dat$masks$test == 1)
dat$Y <- as.matrix(dat$depart)
# dat$Y[is.na(dat$Y)] <- 0
##############################################################################

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
#---------------------------------------------------------------------------------
test_error <<- error_metric$rmse
data_dir <- "./BIXI/results/"

metrics = c(
 "time",
 "lambda.M",
 "lambda.beta",
 "error.test",
 "corr.test",
 "error.train",
 "rank_M",
 "rank_beta",
 "sparse_prop",
 "Prop_explained_xbeta",
 "Prop_explained_M",
 "Prop_unexplained"
)

model_functions = list(
 SImpute_Bixi_Wrapper,
 Mao_Bixi_Wrapper,
 CASMC_0_Bixi_Wrapper,
 #CASMC_1_Bixi_Wrapper,
 CASMC_2_Bixi_Wrapper,
 CASMC_3a_Bixi_Wrapper,
 #CASMC_3b_Bixi_Wrapper,
 Naive_Bixi_Wrapper
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

stopifnot(length(models) == length(model_functions))
perf_means <- perf_stdev <-
 matrix(0,
        length(models),
        length(metrics),
        dimnames = list(models, metrics))
seed = 0

cov_colnames <- c()
statistics <-
 c("Min",
   "1st_Qu",
   "Median",
   "Mean",
   "3rd_Qu",
   "Max",
   "prop_non_zero")
for (var in colnames(dat$X))
 cov_colnames <- c(cov_colnames, paste(var, statistics, sep = "::"))
cov_means <- cov_stdev <-
 matrix(0,
        length(models),
        length(cov_colnames),
        dimnames = list(models, cov_colnames))




for (n in 1:num_replications) {
 seed = seed + n
 set.seed(seed)
 # load data
 if(spatial){
   
 dat <-
  load_bixi_dat(transpose = F, scale_response = T, scale_covariates = F,
                testp = 0.2, validp = 0.2, seed=seed)$model
 dat$X <- dat$spatial_simple
 dat$X <- dat$X |> #[, c(1, 2, 5)] |>
  scalers("minmax")
 }else{
 dat <-
  load_bixi_dat(
   transpose = T,
   scale_response = T,
   scale_covariates = F,
   testp = 0.2,
   validp = 0.2,
   seed = seed
  )$model
 dat$X <- dat$X[, c(1, 2, 5)] |>
   scalers("minmax")
 
 }
 
 
 
 
 
 dat$masks$tr_val = (dat$masks$obs == 1) & (dat$masks$test == 1)
 dat$Y <- as.matrix(dat$depart)
 #-------------------------------------------------------------
 print(n)
 
 for (i in 1:length(model_functions)) {
   print(paste("starting Model ",i, ": ", models[i]))
  results =
   model_functions[[i]](
    dat = dat,
    n_folds = 5,
    lambda.1_grid = seq(1, 0, length = 20),
    lambda.2_grid = seq(.9, 0, length = 20),
    alpha_grid = c(1),
    ncores = 1,
    max_cores = 20,
    weight_function = Mao_weights$uniform
   )
  
  no_cov = ifelse(i == 1 || i == 7, TRUE, FALSE)
  if (!no_cov)
   cov_summs <- (results$cov_summaries |>
                  pivot_longer(cols = everything()) |>
                  t())[2, ] |>
   as.numeric()
  
  results$cov_summaries <- NULL
  results$model = NULL
  results <- results |> as.numeric()
  
  
  new_mean <- update_mean(perf_means[i,], results, n)
  perf_stdev[i,] <-
   update_sse(perf_means[i,], new_mean, perf_stdev[i,],
              results, n)
  perf_means[i,] <- new_mean
  #-----
  # update cov
  if (!no_cov) {
   new_mean <- update_mean(cov_means[i,], cov_summs, n)
   cov_stdev[i,] <-
    update_sse(cov_means[i,], new_mean, cov_stdev[i,],
               cov_summs, n)
   cov_means[i,] <- new_mean
  }
  
 }
 print(perf_means)
 #print(cov_means)
 # print(perf_stdev)
 
 #-------------------------------------
}
for (i in 1:length(model_functions)) {
 perf_stdev[i,] <- get_sd_from_sse(perf_stdev[i,], n)
 cov_stdev[i,] <- get_sd_from_sse(cov_stdev[i,], n)
}
#---------------------
# we now combine them in a dataframe
df_means <- as.data.frame(perf_means)
df_stdev <- as.data.frame(perf_stdev)

df_cov_means <-
 as.data.frame(cov_means) |>
 mutate(Model = rownames(cov_means)) |>
 pivot_longer(
  -Model,
  names_to = c("Variable", "statistic"),
  names_sep = "::",
  values_to = "value"
 ) |>
 pivot_wider(names_from = statistic, values_from = value)


df_cov_stdev <-
 as.data.frame(cov_stdev) |>
 mutate(Model = rownames(cov_stdev)) |>
 pivot_longer(
  -Model,
  names_to = c("Variable", "statistic"),
  names_sep = "::",
  values_to = "value"
 ) |>
 pivot_wider(names_from = statistic, values_from = value)



# Rename columns
colnames(df_means) <-
 paste(colnames(df_means), "mean", sep = "_")
colnames(df_stdev) <- paste(colnames(df_stdev), "sd", sep = "_")
results <- cbind(df_means, df_stdev)
results$B = n
results$model = models

#print(results)
print("Exiting Loop ...")
keyword = ifelse(spatial, "station_", "")


filename = paste0(keyword, "Model_results_",
                  "",
                  ".csv")

write.csv(results,
          file = paste0(data_dir, filename),
          row.names = FALSE)

write.csv(df_cov_means,
          file = paste0(data_dir, paste0(keyword,"cov_means.csv")),
          row.names = FALSE)

write.csv(df_cov_stdev,
          file = paste0(data_dir, paste0(keyword,"cov_std.csv")),
          row.names = FALSE)


print(paste("Results saved to", paste0(data_dir, filename)))
#----------------------------