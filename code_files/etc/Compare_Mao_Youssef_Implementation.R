setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/")
source("./Mao/SMC_functions6.R")
source("./code_files/Mao_import_lib.R")
data_dir = "./Mao/saved_data/"


missingness = 0.9
dim = seq(200,800,200)
results <- data.frame(Dim=paste0(dim,"x",dim), true_rank=NA, 
                      new.alpha=NA, new.lambda.1=NA, new.lambda.2=NA, new.error.test=NA,
                      new.error.all=NA, new.error.B=NA, new.error.beta=NA, new.rank=NA,
                      Ysf.alpha=NA, Ysf.lambda.1=NA, Ysf.lambda.2=NA, Ysf.error.test=NA,
                      Ysf.error.all=NA, Ysf.error.B=NA, Ysf.error.beta=NA, Ysf.rank=NA)

for(i in 1:length(dim)){
   
   #gen.dat <- generate_simulation_data_mao(n1=dim[i],n2=dim[i],m=5,r=10, seed=2023)
   gen.dat <- generate_simulation_data_ysf(2,dim[i],dim[i],5,10, missing_prob = missingness,coll=TRUE)
   results$true_rank[i] = gen.dat$rank
   # fit 1.
   cv.out <- Mao.cv(gen.dat$A, gen.dat$X, gen.dat$Y, gen.dat$W,
                    n_folds=5, 
                    lambda.1_grid = seq(0,5,length=40),
                    lambda.2_grid = seq(.9, 0, length=20),
                    alpha_grid = c(1),#seq(0.992, 1, length=10),
                    numCores = 5,n1n2_optimized = TRUE,theta_estimator = theta_default)
   mao.out <- Mao.fit(gen.dat$Y, gen.dat$X, gen.dat$W, cv.out$best_parameters$lambda.1, 
                      cv.out$best_parameters$lambda.2, cv.out$best_parameters$alpha, 
                      theta_estimator = theta_default)
   # results out of fit 1
   results$new.alpha[i] = cv.out$best_parameters$alpha
   results$new.lambda.1[i] = cv.out$best_parameters$lambda.1
   results$new.lambda.2[i] = cv.out$best_parameters$lambda.2
   results$new.error.test[i] = test_error(mao.out$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$new.error.all[i] = test_error(mao.out$A_hat, gen.dat$A)
   results$new.error.B[i] = test_error(mao.out$B_hat, gen.dat$B)
   results$new.error.beta[i] = test_error(mao.out$beta_hat, gen.dat$beta.x)
   results$new.rank[i] = mao.out$rank
   print(i)
   #---------------
   # Youssef's model
   cv_ratio_SMC_x = SMCfit_cv(gen.dat$A, gen.dat$X, gen.dat$W, gen.dat$Y, nfolds = 5, 
                              tau1_grid = seq(0, 5, length = 40), tau2_grid = seq(0.9, 0, length = 20),
                              alpha_grid = c(1),seed=2023)
   fit_mao_x <- SMCfit(gen.dat$Y, gen.dat$X, cv_ratio_SMC_x$tau_beta_ratio, cv_ratio_SMC_x$tau_svd_ratio, 
                       cv_ratio_SMC_x$alpha_ratio)
   # results out of fit 2
   results$Ysf.alpha[i] = cv_ratio_SMC_x$alpha_ratio
   results$Ysf.lambda.1[i] = cv_ratio_SMC_x$tau_beta_ratio
   results$Ysf.lambda.2[i] = cv_ratio_SMC_x$tau_svd_ratio
   results$Ysf.error.test[i] = test_error(fit_mao_x$Ahat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$Ysf.error.all[i] = test_error(fit_mao_x$Ahat, gen.dat$A)
   results$Ysf.error.B[i] = test_error(fit_mao_x$Bhat, gen.dat$B)
   results$Ysf.error.beta[i] = test_error(fit_mao_x$betahat, gen.dat$beta.x)
   results$Ysf.rank[i] = fit_mao_x$rank
   
   print(results[i,])
   
}

write.csv(results, file = paste0(data_dir, "compare_Mao_Youssef_Implementation_YsfSim_90.csv"),
          row.names = FALSE)
   
   
   