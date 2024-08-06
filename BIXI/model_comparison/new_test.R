
source("./code_files/import_lib.R")
source("./BIXI/data/transform_data_into_mat_other.R")

dat <- load_model_bixi_dat2(time_cov = FALSE,2023,.2, .4, .2)

CASMC_Ridge_Sim_Wrapper(dat)
CASMC_Lasso_Sim_Wrapper(dat)


CASMC_Nuclear_Sim_Wrapper(dat, trace=T)


dat$X


ridge_hpar <- CASMC_Ridge_hparams

out <- CASMC_Ridge_cv(dat$fit_data$train, 
                      dat$X, dat$fit_data$valid, 
                      dat$fit_data$W_valid,
                      dat$fit_data$Y, ridge_hpar,
                      utils$error_metric$rmse, max_cores = 20,
                      trace = T, track = T)

nuclear_hpar <- CASMC_Nuclear_hparams

out <- out_nuc <- CASMC_Nuclear_cv(dat$fit_data$train, dat$X, dat$fit_data$valid, dat$fit_data$W_valid,
                                   dat$fit_data$Y, nuclear_hpar,
                                   trace = T, track = T)

out <- readRDS("./saved_data/movie_lens_nuclear_fit.rds")
# saveRDS(out, "./saved_data/movie_lens_nuclear_fit.rds")


out$hparams
fiti <- out
fit. = fiti$fit
# get estimates and validate
fit.$M = utils$unsvd(fit.)
fit.$beta = utils$unsvd(fit.$beta)
fit.$estimates = fit.$M + dat$X %*% fit.$beta

apply(fit.$beta, 1, summary) %>% as.data.frame() %>% t() %>%  round(3)


results = list(model = "CASMC-Nuclear")
results$lambda.beta = fiti$hparams$lambda.beta
results$lambda.M = fiti$hparams$lambda.M
prepare_output(Sys.time(), fit.$estimates,
               dat$O, dat$W, NA, fit.$beta, NA, fit.$M, NULL)  


lasso_hpar <- CASMC_Lasso_hparams
lasso_hpar$beta$lambda.max <- 60


out <- CASMC_Lasso_cv(dat$fit_data$train, dat$X, dat$fit_data$valid,
                      dat$fit_data$W_valid,
                      dat$fit_data$Y, lasso_hpar, max_cores = 20,
                      trace = T)

CASMC_Ridge_Sim_Wrapper(dat)
CASMC_Lasso_Sim_Wrapper(dat)
CASMC_Nuclear_Sim_Wrapper(dat)





