library(tidyverse)
library(dplyr)
library(magrittr)
library(microbenchmark)


source("./code_files/import_lib.R")

#---------------------------------

find_lasso_max_param(
  y_train = dat$fit_data$train,
  X = dat$fit_data$Xq,
  y_valid = NULL,#dat$fit_data$valid,
  W_valid = NULL,#dat$fit_data$W_valid,
  y = NULL,#dat$fit_data$Y,
  maxit = 100,
  verbose=5)
#---------------------------------------
dat <- generate_simulated_data(600, 700, 10, 12, .8, F,
                               informative_cov_prop = .3,
                               prepare_for_fitting = T)
hpar <- CAMC_Lasso_hparams
hpar$beta$n.lambda = 80
hpar$beta$lambda_max = 4
fit5 <- CAMC_Lasso_cv(
  y_train = dat$fit_data$train,
  X = dat$fit_data$Xq,
  y_valid = dat$fit_data$valid,
  W_valid = dat$fit_data$W_valid,
  y = dat$fit_data$Y,
  hpar = hpar,
  verbose = 1,
  max_cores = 6
)

hpar$M

Mao_Sim_Wrapper(dat)$time
SImpute_Sim_Wrapper(dat, hpar)$time
CAMC_Sim_Wrapper(dat, hpar=hpar, max_cores = 10, verbose=0)$time
Naive_Sim_Wrapper(dat)$time
#-----------------------------------

#---------------------------------

#-------------------------------------------

microbenchmark(
  
  meth1 = {
  },
  meth2 = {
      },
  times = 100L
  
)