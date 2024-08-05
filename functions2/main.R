
source("./code_files/import_lib.R")

dat <- generate_simulated_data(300, 350, r= 10, k = 6, missing_prob = 0.8,  informative_cov_prop = .5,
                               mar_sparse = T,
                               prepare_for_fitting = T)

ridge_hpar <- CASMC_Ridge_hparams
out <- CASMC_Ridge_cv(dat$fit_data$train, dat$X, dat$fit_data$valid, dat$fit_data$W_valid,
                      dat$fit_data$Y, ridge_hpar, utils$error_metric$rmse, max_cores = 20,
                      trace = T, track = T)

nuclear_hpar <- CASMC_Nuclear_hparams

out <- CASMC_Nuclear_cv(dat$fit_data$train, dat$X, dat$fit_data$valid, dat$fit_data$W_valid,
                      dat$fit_data$Y, nuclear_hpar,
                      trace = T, track = T)


lasso_hpar <- CASMC_Lasso_hparams
lasso_hpar$beta$lambda.max <- 60


out <- CASMC_Lasso_cv(dat$fit_data$train, dat$X, dat$fit_data$valid,
                      dat$fit_data$W_valid,
                      dat$fit_data$Y, lasso_hpar, max_cores = 20,
                      trace = T)

CASMC_Ridge_Sim_Wrapper(dat)
CASMC_Lasso_Sim_Wrapper(dat)
CASMC_Nuclear_Sim_Wrapper(dat)
