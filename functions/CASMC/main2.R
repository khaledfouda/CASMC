# here we compare/test our implementation
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

# 0. prepare the data

dat <-  generate_simulation_rows(800,900, r = 10, k = 10, missing_prob = 0.9, coll=FALSE,
                                 prepare_for_fitting = TRUE,
                                  half_discrete = FALSE, informative_cov_prop = 1, seed = 2023)
X_r = reduced_hat_decomp(dat$X, 1e-2)


best_fit2 = CASMC_cv_rank(Y_train, X_r, Y_valid, W_valid,  y=y, 
                                      lambda.beta.grid =  "default",#"default",#seq(0,5,length.out=20), 
                                      max_cores = 30,
                                      trace=F, thresh=1e-6,n.lambda = 30, rank.limit = 20, track_beta = T) 
