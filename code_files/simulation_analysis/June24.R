dat <-
 generate_simulation_rows(
  600,
  700,
  r = 10,
  k = 10, 
  missing_prob = 0.9,
  coll = T,
  prepare_for_fitting = TRUE,
  half_discrete = FALSE,
  informative_cov_prop = 1,mar_sparse = T,
  mv_beta = T,
  seed = 2023
 )


fiti <- CASMC2_cv(
 y_train = dat$fit_data$train,
 X = dat$X,
 y_valid = dat$fit_data$valid,
 W_valid = dat$fit_data$W_valid,
 y = dat$fit_data$Y,
 error_function = error_metric$rmse,
 warm = NULL,
 quiet = F,
 trace = F,
 track = T,
 rank.beta.init = 1,
 rank.beta.step = 2,
 lambda.beta.grid = "default1",
 max_cores = max_cores,
 seed = NULL,
)

fiti


system.time({for(i in 1:5000) x=t( dd$d * t(dd$u)) })
