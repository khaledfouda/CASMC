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
          showres=T) {
  error_function_name <-
   names(which(sapply(error_metric, identical, test_error)))
  
  if (mao) {
   M = outp$M
   O = outp$estimates
  } else{
   M = outp$u %*% (outp$d * t(outp$v))
   O = M + inp$X %*% outp$beta
  }
  error_beta <- round(test_error(outp$beta, inp$beta), 5)
  error_M <- round(test_error(M, inp$M), 5)
  error_test <-
   round(test_error(O[inp$W == 0], inp$O[inp$W == 0]), 5)
  error_train <-
   round(test_error(O[inp$W != 0], inp$O[inp$W != 0]), 5)
  
  rank_beta <- qr(outp$beta)$rank
  rank_M <- qr(M)$rank
  
  result_df <- data.frame(
   model = name,
   Metric = c(paste0(error_function_name, "(rank)")),
   Beta = sprintf("%.5f(%d)", error_beta, rank_beta),
   M = sprintf("%.5f(%d)", error_M, rank_M),
   Test = sprintf("%.5f", error_test),
   Train = sprintf("%.5f", error_train)
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

metric = error_metric$mae
kable( arrange(rbind(
print_performance(dat, fit_rank$fit, metric, F, "CASMC(Rank)",F),
print_performance(dat, fit_l2$fit, metric, F, "CASMC(L2)",F),
print_performance(dat, mao.fit$fit, metric, TRUE, "Mao(Uni)",F),
print_performance(dat, mao.fitB$fit, metric, TRUE, "Mao(Binom)",F)
),Test), format = "simple")

mao.fitB <- Mao.cv(
 Y = dat$Y,
 X = dat$X,
 W = dat$W,
 n_folds = 5,
 lambda.1_grid = seq(0, 1, length = 20),
 lambda.2_grid = seq(0.9, 0.1, length = 20),
 alpha_grid = seq(0.992, 1, length = 20),
 seed = 2023,
 numCores = 1,
 n1n2_optimized = TRUE,
 test_error = error_metric$rmse,
 theta_estimator = Mao_weights$binomial,
 sequential = FALSE
)
mao.fitB$best_parameters

rbind(
dat$beta[,1] |> round(4),
mao.fitB$fit$beta[,1] |> round(4),
fit_rank$fit$beta[,1] |> round(4)
)
