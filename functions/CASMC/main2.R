# here we compare/test our implementation
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

# 0. prepare the data

print_performance <-
 function(inp,
          outp,
          test_error = error_metric$rmse,
          mao = FALSE,
          name = "",
          showres=T,
          rdig = 4) {
  error_function_name <-
   names(which(sapply(error_metric, identical, test_error)))
  
  if (mao) {
   M = outp$M
   O = outp$estimates
  } else{
   M = outp$u %*% (outp$d * t(outp$v))
   O = M + inp$X %*% outp$beta
  }
  error_beta <- round(test_error(outp$beta, inp$beta), rdig)
  error_M <- round(test_error(M, inp$M), rdig)
  error_test <-
   round(test_error(O[inp$W == 0], inp$O[inp$W == 0]), rdig)
  error_train <-
   round(test_error(O[inp$W != 0], inp$O[inp$W != 0]), rdig)
  error_Y <- 
   round(test_error(O[inp$W != 0], inp$Y[inp$W != 0]), rdig)
  
  rank_beta <- qr(outp$beta)$rank
  rank_M <- qr(M)$rank
  
  result_df <- data.frame(
   Metric = c(paste0(error_function_name, "(rank)")),
   model = name,
   Beta = sprintf(paste0("%.",rdig,"f(%2d)"), error_beta, rank_beta),
   M = sprintf(paste0("%.",rdig,"f(%2d)"), error_M, rank_M),
   test = sprintf(paste0("%.",rdig,"f"), error_test),
   train = sprintf(paste0("%.",rdig,"f"), error_train),
   Y = sprintf(paste0("%.",rdig,"f"), error_Y)
  )
  
  if(showres){
  print(knitr::kable(result_df, format = "simple" ))
  }else
 return(result_df)
 }


dat <-
 generate_simulation_rows(
  800,
  900,
  r = 10,
  k = 10,
  missing_prob = 0.9,
  coll = FALSE,
  prepare_for_fitting = TRUE,
  half_discrete = FALSE,
  informative_cov_prop = 0.7,
  seed = 2023
 )
X_r = reduced_hat_decomp(dat$X, 1e-2)

fit_rank <- CASMC_cv_rank(
 y_train = dat$fit_data$train,
 X = dat$X,
 y_valid = dat$fit_data$valid,
 W_valid = dat$fit_data$W_valid,
 y = dat$fit_data$Y,
 error_function = error_metric$rmse,
 lambda.factor = 1 / 4,
 lambda.init = NULL,
 n.lambda = 20,
 rank.init = 2,
 rank.limit = 30,
 rank.step = 2,
 pct = 0.98,
 lambda.a = 0,
 S.a = NULL,
 lambda.b = 0,
 S.b = NULL,
 early.stopping = 1,
 thresh = 1e-6,
 maxit = 100,
 trace = FALSE,
 print.best = TRUE,
 quiet = FALSE,
 warm = NULL,
 rank_x = X_r$rank,
 r_min = 0,
 r_max = X_r$rank,
 track_r = TRUE,
 max_cores = 20,
 seed = 2023
)





fit_l2 <- CASMC_cv_L2(
 y_train = dat$fit_data$train,
 X = dat$X,
 y_valid = dat$fit_data$valid,
 W_valid = dat$fit_data$W_valid,
 y = dat$fit_data$Y,
 error_function = error_metric$rmse,
 lambda.factor = 1 / 4,
 lambda.init = NULL,
 n.lambda = 20,
 rank.init = 2,
 rank.limit = 30,
 rank.step = 2,
 pct = 0.98,
 lambda.a = 0,
 S.a = NULL,
 lambda.b = 0,
 S.b = NULL,
 early.stopping = 1,
 thresh = 1e-6,
 maxit = 100,
 trace = FALSE,
 print.best = TRUE,
 quiet = FALSE,
 warm = NULL,
 lambda.beta.grid = "default",
 track_beta = TRUE,
 max_cores = 23,
 seed = 2023
)
fit_l2_ <- CASMC_cv_L2(
 y_train = dat$fit_data$train,
 X = dat$X,
 y_valid = dat$fit_data$valid,
 W_valid = dat$fit_data$W_valid,
 y = dat$fit_data$Y,
 error_function = error_metric$rmse,
 lambda.factor = 1 / 4,
 lambda.init = NULL,
 n.lambda = 20,
 rank.init = 2,
 rank.limit = 30,
 rank.step = 2,
 pct = 0.98,
 lambda.a = 0,
 S.a = NULL,
 lambda.b = 0,
 S.b = NULL,
 early.stopping = 1,
 thresh = 1e-6,
 maxit = 2,
 trace = FALSE,
 print.best = TRUE,
 quiet = FALSE,
 warm = NULL,
 lambda.beta.grid = "default",
 track_beta = TRUE,
 max_cores = 23,
 seed = 2023
)

mao.fitB <- Mao.cv(
 Y = dat$Y,
 X = dat$X,
 W = dat$W,
 n_folds = 5,
 lambda.1_grid = seq(0, 1, length = 20),
 lambda.2_grid = seq(0.9, 0.1, length = 20),
 alpha_grid = c(1),#seq(0.992, 1, length = 5),
 seed = 2023,
 numCores = 1,
 n1n2_optimized = TRUE,
 test_error = error_metric$rmse,
 theta_estimator = Mao_weights$binomial,
 sequential = FALSE
)




mao.fit <- Mao.cv(
 Y = dat$Y,
 X = dat$X,
 W = dat$W,
 n_folds = 5,
 lambda.1_grid = seq(0, 1, length = 20),
 lambda.2_grid = seq(0.9, 0.1, length = 20),
 alpha_grid = c(1),#seq(0.992, 1, length = 5),
 seed = 2023,
 numCores = 1,
 n1n2_optimized = TRUE,
 test_error = error_metric$rmse,
 theta_estimator = Mao_weights$uniform,
 sequential = FALSE
)

mao.fitB$best_parameters

naive_model <- naive_fit(dat$Y, dat$X)


metric = error_metric$rmse
results <- rbind(
 print_performance(dat, fit_rank$fit, metric, F, "CASMC(Rank)",F,3),
 print_performance(dat, fit_l2$fit, metric, F, "CASMC(L2)",F,3),
 print_performance(dat, mao.fit$fit, metric, TRUE, "Mao(Uni)",F,3),
 print_performance(dat, mao.fitB$fit, metric, TRUE, "Mao(Binom)",F,3),
 #print_performance(dat, fit_rank_$fit, metric, F, "CASMC(Rank) 1 iter",F,3),
 print_performance(dat, fit_l2_$fit, metric, F, "CASMC(L2) 1 iter",F, 3),
 print_performance(dat, naive_model, metric, TRUE, "Naive", F, 3))

results  |> 
 arrange(Beta) |> 
 mutate(Beta = paste0("[",1:nrow(results),"]",Beta)) |> 
 arrange(M) |> 
 mutate(M = paste0("[",1:nrow(results),"]",M)) |> 
 arrange(Y) |> 
 mutate(Y = paste0("[",1:nrow(results),"]",Y)) |> 
 select(-Y) |> 
 arrange(train) |> 
 mutate(train = paste0("[",1:nrow(results),"]",train)) |> 
 arrange(test) |> 
 mutate(test = paste0("[",1:nrow(results),"]",test))  ->
 results


kable( results, format = "simple")
#----------------------------------------------------------------------------------



fit_l2$fit$M <- unsvd(fit_l2$fit)
fit_rank$fit$M <- unsvd(fit_rank$fit)


abs((dat$X %*% dat$beta)[1,1:6] -(dat$X %*%  mao.fit$fit$beta)[1,1:6]) |> round(2) - 
abs((dat$X %*% dat$beta)[1,1:6] - (dat$X %*% fit_l2$fit$beta)[1,1:6]) |> round(2)


abs(dat$M[1,1:6] - mao.fit$fit$M[1,1:6]) |> round(2) - 
 abs(dat$M[1,1:6] - fit_l2$fit$M[1,1:6]) |> round(2)


 abs((dat$Y)[1,1:6] -(dat$X %*%  mao.fit$fit$beta+ mao.fit$fit$M)[1,1:6]) |> round(2) -
 abs((dat$Y)[1,1:6] - (dat$X %*% fit_l2$fit$beta+fit_l2$fit$M)[1,1:6]) |> round(2)





