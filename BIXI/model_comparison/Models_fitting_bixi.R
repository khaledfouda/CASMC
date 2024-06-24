


setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
source("./code_files/import_lib.R")
source("./BIXI/data-raw/bixi_data.R")

dat <-
 load_bixi_dat(transpose = T, scale_response = T, scale_covariates = F,
               testp = 0.2, validp = 0.2, seed=2023)$model
dat$X <- dat$X[,c(1,2,5)] |> 
 scalers("minmax")

dat$masks$tr_val = (dat$masks$obs == 1) & (dat$masks$test == 1)

