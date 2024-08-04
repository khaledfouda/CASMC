utils <- new.env()

source("./code_files/import_lib.R")

dat <- generate_simulated_data(300, 350, r= 10, k = 6, missing_prob = 0.8, 
                               prepare_for_fitting = T)


out <- CASMC_Ridge_cv(dat$fit_data$train, dat$X, dat$fit_data$valid, dat$fit_data$W_valid,
                      dat$fit_data$Y, utils$error_metric$rmse, max_cores = 20)


out <- CASMC_Nuclear_cv(dat$fit_data$train, dat$X, dat$fit_data$valid, dat$fit_data$W_valid,
                      dat$fit_data$Y, beta_cv_param = list(
                              rank.init = 2,
                              rank.limit = qr(dat$X)$rank,
                              rank.step = 2,
                              pct = 0.98,
                              lambda.factor = 1,
                              lambda.init = NULL,
                              n.lambda = 20,
                              early.stopping = 1
                      ),
                      trace = T, track = T)
