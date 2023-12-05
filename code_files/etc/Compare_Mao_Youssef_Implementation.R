

missingness = 0.9

compare_and_save <- function(missingness,coll=TRUE, 
                             lambda.1_grid = seq(0,3,length=20),
                             lambda.2_grid = seq(.9, 0, length=20),
                             alpha_grid = seq(0.992, 1, length=10), plot=FALSE, tofile=FALSE, graph_label="",ncores=2){
   
   setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/")
   source("./Mao/SMC_functions6.R")
   source("./code_files/Mao_import_lib.R")
   data_dir = "./Mao/saved_data/"
   stopifnot(missingness %in% c(0,0.8, 0.9))
   dim = seq(400,1000,200)
   results <- data.frame(Dim=paste0(dim,"x",dim), true_rank=NA, 
                         new.alpha=NA, new.lambda.1=NA, new.lambda.2=NA, new.error.test=NA,
                         new.error.all=NA, new.error.B=NA, new.error.beta=NA, new.rank=NA, new.time=NA,
                         Ysf.alpha=NA, Ysf.lambda.1=NA, Ysf.lambda.2=NA, Ysf.error.test=NA,
                         Ysf.error.all=NA, Ysf.error.B=NA, Ysf.error.beta=NA, Ysf.rank=NA, Ysf.time=NA)
   
   for(i in 1:length(dim)){
      
      if(missingness == 0){
         gen.dat <- generate_simulation_data_mao(n1=dim[i],n2=dim[i],m=5,r=10, seed=2023)
      }else
         gen.dat <- generate_simulation_data_ysf(2,dim[i],dim[i],10,10, missing_prob = missingness,coll=coll)
      results$true_rank[i] = gen.dat$rank
      set.seed(2023)
      # fit 1.
      start_time = Sys.time()
      cv.out <- Mao.cv(gen.dat$A, gen.dat$X, gen.dat$Y, gen.dat$W,
                       n_folds=5, 
                       lambda.1_grid = lambda.1_grid,
                       lambda.2_grid = lambda.2_grid,
                       alpha_grid = alpha_grid,
                       numCores = ncores,n1n2_optimized = TRUE,theta_estimator = theta_default)
      mao.out <- Mao.fit(gen.dat$Y, gen.dat$X, gen.dat$W, cv.out$best_parameters$lambda.1, 
                         cv.out$best_parameters$lambda.2, cv.out$best_parameters$alpha, 
                         theta_estimator = theta_default)
      results$new.time[i] = round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
      # results out of fit 1
      results$new.alpha[i] = cv.out$best_parameters$alpha
      results$new.lambda.1[i] = cv.out$best_parameters$lambda.1
      results$new.lambda.2[i] = cv.out$best_parameters$lambda.2
      results$new.error.test[i] = test_error(mao.out$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
      results$new.error.all[i] = test_error(mao.out$A_hat, gen.dat$A)
      results$new.error.B[i] = test_error(mao.out$B_hat, gen.dat$B)
      results$new.error.beta[i] = test_error(mao.out$beta_hat, gen.dat$beta)
      results$new.rank[i] = mao.out$rank
      print(i)
      #---------------
      # Youssef's model
      set.seed(2023)
      start_time = Sys.time()
      cv_ratio_SMC_x = SMCfit_cv(gen.dat$A, gen.dat$X, gen.dat$W, gen.dat$Y, nfolds = 5, 
                                 tau1_grid = lambda.1_grid, tau2_grid = lambda.2_grid,
                                 alpha_grid = alpha_grid,seed=2023)
      fit_mao_x <- SMCfit(gen.dat$Y, gen.dat$X, cv_ratio_SMC_x$tau_beta_ratio, cv_ratio_SMC_x$tau_svd_ratio, 
                          cv_ratio_SMC_x$alpha_ratio)
      # results out of fit 2
      results$Ysf.time[i] =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
      results$Ysf.alpha[i] = cv_ratio_SMC_x$alpha_ratio
      results$Ysf.lambda.1[i] = cv_ratio_SMC_x$tau_beta_ratio
      results$Ysf.lambda.2[i] = cv_ratio_SMC_x$tau_svd_ratio
      results$Ysf.error.test[i] = test_error(fit_mao_x$Ahat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
      results$Ysf.error.all[i] = test_error(fit_mao_x$Ahat, gen.dat$A)
      results$Ysf.error.B[i] = test_error(fit_mao_x$Bhat, gen.dat$B)
      results$Ysf.error.beta[i] = test_error(fit_mao_x$betahat, gen.dat$beta)
      results$Ysf.rank[i] = fit_mao_x$rank
      #-------------------------------------
      # soft Impute model
      
      
     
      
      
      # saving plots to disk
      if(plot==TRUE){
         filename = paste0(graph_label,"_theta" ,missingness, c("_A","_beta", "_B"), "_dim",dim[i],"_coll",coll)
         plot_actual.vs.estimated.v2(gen.dat$A[gen.dat$W==0], mao.out$A_hat[gen.dat$W==0], fit_mao_x$Ahat[gen.dat$W==0],
                                     "New", "Old", "A", expression(hat("A")), "",sample.size = 10000, tofile=tofile, filename=filename[1])
         plot_actual.vs.estimated.v2(gen.dat$beta, mao.out$beta_hat, fit_mao_x$betahat,
                                     "New", "Old", "Beta", expression(hat("Beta")), "",sample.size = 10000, tofile=tofile, filename=filename[2])
         plot_actual.vs.estimated.v2(gen.dat$B, mao.out$B_hat, fit_mao_x$Bhat,
                                     "New", "Old", "B", expression(hat("B")), "",sample.size = 10000, tofile=tofile, filename=filename[3])
      }
      print(results[i,])
      
   }
   if(missingness == 0){
      filename = "compare_Mao_Youssef_Implementation_MaoSim.csv"
   }else
      filename = paste0("compare_Mao_Youssef_Implementation_YsfSim_",round(missingness*100),
                        "_coll_",coll, ".csv")
   
   write.csv(results, file = paste0(data_dir, filename), row.names = FALSE)
   
   #----------------------------
}

alpha_grid = c(1)
ncores = 1

compare_and_save(0.8, TRUE, plot = TRUE, tofile = TRUE,
                 lambda.1_grid = seq(0,2,length=20),lambda.2_grid = seq(.9, 0, length=20),
                 alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0.9, TRUE, plot = TRUE, tofile = TRUE,
                 lambda.1_grid = seq(0,2,length=20),lambda.2_grid = seq(.9, 0, length=20),
                 alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0.8, FALSE, plot = TRUE, tofile = TRUE,
                 lambda.1_grid = seq(0,2,length=20),lambda.2_grid = seq(.9, 0, length=20),
                 alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0.9, FALSE, plot = TRUE, tofile = TRUE,
                 lambda.1_grid = seq(0,2,length=20),lambda.2_grid = seq(.9, 0, length=20),
                 alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0,  FALSE, plot = TRUE, tofile = TRUE,
                 lambda.1_grid = seq(0,2,length=20), lambda.2_grid = seq(.9, 0, length=20), 
                 alpha_grid = alpha_grid, ncores=ncores)

#------------------------------------------
image_files <- function(miss,coll, dim, var) 
   paste0("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/Mao/saved_data/rds_plots/_theta",
          miss,"_",var,"_dim",dim,"_coll",coll,".rds")
