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
# working with bixi data
source("./BIXI/data/transform_data_into_mat.R")

dat <- load_model_bixi_dat(time_cov = FALSE)
hpar <- CAMC_Lasso_hparams
hpar$beta$n.lambda = 80
CAMC_Bixi_Wrapper(dat, hpar=hpar, train_on_all = TRUE)
Naive_Bixi_Wrapper(dat)
SImpute_Bixi_Wrapper(dat, hpar)
Mao_Bixi_Wrapper(dat)



fit5 <- CAMC_Lasso_cv(
  y_train = dat$splits$train,
  X = dat$splits$Xq,
  y_valid = dat$splits$valid@x,
  W_valid = dat$masks$valid,
  y = dat$splits$Y,
  hpar = hpar,
  verbose = 1,
  max_cores = 6
)

dat$X[1:5,]
dat$Y[1:5,1:5]
dat$splits$Y[1:5,1:5]


fit$beta %>%
  t() %>%
  as.data.frame() %>%
  round(10) %>%
  mutate(a = rowSums(across(everything(), ~.x))) %>%
  filter(a!=0) %>%
  round(4)

sum(fit$beta==0) / length(fit$beta)

apply(round(fit$beta,5), 1, function(x) sum(x==0) / length(x))
apply(round(fit$Rbeta,5), 1, function(x) sum(x==0) / length(x))


#---------------------------------

#-------------------------------------------

microbenchmark(
  
  meth1 = {
  },
  meth2 = {
      },
  times = 100L
  
)