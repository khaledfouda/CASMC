require(foreach)
require(doParallel)

compare_Mao <- function(gen.dat, lambda.1_grid, lambda.2_grid, alpha_grid, numCores=1, n_folds=5){
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
   
   results = list(model = "Mao")
   results$time = round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
   results$alpha = cv.out$best_parameters$alpha
   results$lambda.1 = cv.out$best_parameters$lambda.1
   results$lambda.2 = cv.out$best_parameters$lambda.2
   results$error.test = test_error(mao.out$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$error.all = test_error(mao.out$A_hat, gen.dat$A)
   results$error.B = test_error(mao.out$B_hat, gen.dat$B)
   results$error.beta = test_error(mao.out$beta_hat, gen.dat$beta)
   results$rank = mao.out$rank
   results
}

compare_softImpute_orig <- function(gen.dat, valid.dat){
   start_time = Sys.time()
   sout <- simpute.orig(valid.dat$Y_train, valid.dat$W_valid,
                        gen.dat$Y, trace=FALSE, rank.limit = 30,print.best=FALSE, rank.step = 4)
   results = list(model = "SoftImpute_Orig")
   results$time =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
   results$alpha = NA
   results$lambda.1 = NA
   results$lambda.2 = sout$lambda
   results$error.test = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$error.all = test_error(sout$A_hat, gen.dat$A)
   results$error.B = NA
   results$error.beta =NA
   results$rank = sout$rank_A 
   results
}

compare_softImpute_cov <- function(gen.dat, valid.dat){
   start_time = Sys.time()
   sout <- simpute.cov.cv(valid.dat$Y_train, gen.dat$X, valid.dat$W_valid, valid.dat$Y_valid,
                          trace=F, rank.limit = 30, 
                          print.best=FALSE, rank.step=4, type="als", lambda1=0, tol=2)
   results = list(model = "SoftImpute_Cov")
   results$time =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
   results$alpha = NA
   results$lambda.1 = 0
   results$lambda.2 = sout$lambda
   results$error.test = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$error.all = test_error(sout$A_hat, gen.dat$A)
   results$error.B = test_error(sout$B_hat, gen.dat$B)
   results$error.beta = test_error(sout$beta_hat, gen.dat$beta)
   results$rank = sout$rank_A
   results
}

compare_softImpute_L2 <- function(gen.dat, valid.dat, lambda.1_grid){
   
   start_time = Sys.time()
   sout <- simpute.cov.cv(valid.dat$Y_train, gen.dat$X, valid.dat$W_valid,
                          valid.dat$Y_valid, trace=FALSE, rank.limit = 30, 
                          print.best=FALSE, rank.step=4, type="als", lambda1=0, tol=2)
   sout <- simpute.cov.cv.lambda1(valid.dat$Y_train, gen.dat$X,valid.dat$W_valid,
                                  valid.dat$Y_valid, sout$lambda, sout$rank.max, print.best = FALSE,
                                  trace=FALSE, lambda1.grid =lambda.1_grid ,n1n2 = 1, warm=NULL)
   
   results = list(model = "SoftImpute_L2")
   results$time =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
   results$alpha = NA
   results$lambda.1 = sout$lambda1
   results$lambda.2 = sout$lambda2
   results$error.test = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$error.all = test_error(sout$A_hat, gen.dat$A)
   results$error.beta = test_error(sout$beta_hat, gen.dat$beta)
   results$error.B = test_error(sout$B_hat, gen.dat$B)
   results$rank = sout$rank_A
   results
}

compare_softImpute_Kfold <- function(gen.dat, lambda.1_grid){
   
   start_time = Sys.time()
   sout <- simpute.cov.kfold(gen.dat$Y, gen.dat$X, gen.dat$W, n_folds = 3, print.best = FALSE,
                             trace=FALSE, rank.limit = 30, lambda1=0,n1n2 = 1, warm=NULL,tol = 2)
   sout <- simpute.cov.kfold.lambda1(gen.dat$Y, gen.dat$X, gen.dat$W, sout$lambda2, n_folds = 3, print.best = FALSE, 
                                     trace=FALSE,lambda1.grid = lambda.1_grid ,n1n2 = 1, warm=NULL,
                                     J=c(sout$J))
   results = list(model = "SoftImpute_Kfold")
   results$time =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
   results$alpha = NA
   results$lambda.1 = sout$lambda1
   results$lambda.2 = sout$lambda2
   results$error.test = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$error.all = test_error(sout$A_hat, gen.dat$A)
   results$error.beta = test_error(sout$beta_hat, gen.dat$beta)
   results$error.B = test_error(sout$B_hat, gen.dat$B)
   results$rank = sout$rank_A
   results
}




compare_softImpute_splr <- function(gen.dat, valid.dat, splr.dat){
   
   start_time = Sys.time()
   
   best_fit = simpute.cov.cv_splr_no_patience(valid.dat$Y_train,
                                              splr.dat$svdH, valid.dat$Y_valid,
                                              valid.dat$W_valid,warm = NULL,
                                              trace = F, rank.limit=30,rank.step=4,patience = 1,
                                              rank.init = 2, lambda.factor = 1/4, n.lambda = 30)
   yfill = gen.dat$Y
   fits = best_fit$best_fit
   best_fit$B_hat = fits$u %*% (fits$d * t(fits$v))
   yfill[gen.dat$Y==0] = (best_fit$B_hat)[gen.dat$Y==0]
   beta_partial = MASS::ginv(t(splr.dat$X) %*% splr.dat$X) %*% t(splr.dat$X)
   fits.out = list(u=fits$u, d=fits$d, v=fits$v, beta.estim=beta_partial %*% yfill)
   
   
   sout <- simpute.als.cov(gen.dat$Y, splr.dat$X, beta_partial,J = best_fit$rank_B, thresh =  1e-6,
                           lambda= best_fit$lambda,
                           trace.it = F,warm.start = fits.out, maxit=100)
   
   sout$B_hat = sout$u %*% (sout$d * t(sout$v))
   sout$A_hat = sout$B_hat  + splr.dat$X %*% sout$beta.estim
   sout$beta_hat = sout$beta.estim
   
   
   results = list(model = "SoftImpute_splr")
   results$time =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
   results$alpha = NA
   results$lambda.1 = NA
   results$lambda.2 = sout$lambda
   results$error.test = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$error.all = test_error(sout$A_hat, gen.dat$A)
   results$error.beta = test_error(sout$beta_hat, gen.dat$beta)
   results$error.B = test_error(sout$B_hat, gen.dat$B)
   results$rank = qr(sout$A_hat)$rank
   results
}

compare_softImpute_splr_Kfold <- function(gen.dat, valid.dat, splr.dat, n_folds=10){
   
   start_time = Sys.time()
   
   best_fit = simpute.cov.Kf_splr_no_patience_v2(gen.dat$Y, splr.dat$svdH, gen.dat$W, n_folds=n_folds,
                                                 trace = F, rank.init=2, lambda.factor = 1/4,n.lambda = 30,
                                                 rank.limit=30,rank.step=4,patience = 1)
   
   
   yfill = gen.dat$Y
   fits = best_fit$best_fit
   best_fit$B_hat = fits$u %*% (fits$d * t(fits$v))
   yfill[gen.dat$Y==0] = (best_fit$B_hat)[gen.dat$Y==0]
   beta_partial = MASS::ginv(t(splr.dat$X) %*% splr.dat$X) %*% t(splr.dat$X)
   fits.out = list(u=fits$u, d=fits$d, v=fits$v, beta.estim=beta_partial %*% yfill)
   
   
   sout <- simpute.als.cov(gen.dat$Y, splr.dat$X, beta_partial,J = best_fit$rank_B, thresh =  1e-6,
                           lambda= best_fit$lambda,
                           trace.it = F,warm.start = fits.out, maxit=100)
   
   sout$B_hat = sout$u %*% (sout$d * t(sout$v))
   sout$A_hat = sout$B_hat  + splr.dat$X %*% sout$beta.estim
   sout$beta_hat = sout$beta.estim
   
   
   results = list(model = "SoftImpute_splr_Kfold")
   results$time =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
   results$alpha = NA
   results$lambda.1 = NA
   results$lambda.2 = sout$lambda
   results$error.test = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$error.all = test_error(sout$A_hat, gen.dat$A)
   results$error.beta = test_error(sout$beta_hat, gen.dat$beta)
   results$error.B = test_error(sout$B_hat, gen.dat$B)
   results$rank = qr(sout$A_hat)$rank
   results
}



compare_and_save <- function(missingness,coll=TRUE,  n_folds=5,
                             dim = seq(400,1000,200), ncovariates=10,
                             lambda.1_grid = seq(0,3,length=20),
                             lambda.2_grid = seq(.9, 0, length=20),
                             alpha_grid = seq(0.992, 1, length=10),ncores=1){
   
   ncores = min(ncores, length(dim))
   if(ncores > 1)
      registerDoParallel(cores=ncores)
   
   print(paste("Running on",ncores,"cores."))
   setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
   source("./code_files/import_lib.R")
   data_dir = "./saved_data/"
   stopifnot(missingness %in% c(0,0.8, 0.9))
   final.results = data.frame()
   
   #for(i in 1:length(dim)){
   
   final.results <- foreach(i = 1:length(dim), .combine='rbind') %dopar% {
    
      set.seed(2023)
      results <- data.frame()
      if(missingness == 0){
         gen.dat <- generate_simulation_data_mao(n1=dim[i],n2=dim[i],m=ncovariates,r=ncovariates, seed=2023)
      }else
         gen.dat <- generate_simulation_data_ysf(2,dim[i],dim[i],ncovariates,ncovariates, 
                                                 missing_prob = missingness,coll=coll)
      #-------------------------------------------------------------------------------------
      # validation set to be used for the next two models
      valid.dat = list()
      valid.dat$W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
      valid.dat$Y_train <- gen.dat$Y * valid.dat$W_valid
      valid.dat$Y_valid <- gen.dat$Y[valid.dat$W_valid==0]
      #---------------------------------------------------------
      # SPLR data
      splr.dat = reduced_hat_decomp(gen.dat$X, 1e-2)
      #----------------------------------------------------------------------
      # fit 1. Mao
      print(i)
      results = rbind(results, compare_Mao(gen.dat, lambda.1_grid, lambda.2_grid,
                                           alpha_grid,1,n_folds))
      print(".")
      #----------------------------------------------------------
      # soft Impute model without covariates
      results = rbind(results, compare_softImpute_orig(gen.dat, valid.dat))
      print("..")
      #----------------------------------------------------------------------------
      # soft Impute model with covariates
      results = rbind(results, compare_softImpute_cov(gen.dat, valid.dat)) 
      print("...")
      #-------------------------------------------------------------------------------------
      # Soft Impute with Covariates and With L2 regularization on the covariates
      results = rbind(results, compare_softImpute_L2(gen.dat, valid.dat, lambda.1_grid))
      print("....")
      #-------------------------------------------------------------------------------------
      # Soft Impute with Covariates and With L2 regularization on the covariates and K-fold cross-validation
      results = rbind(results, compare_softImpute_Kfold(gen.dat, lambda.1_grid))
      print(".....")
      #--------------------------------------------------------------------------------
      results = rbind(results, compare_softImpute_splr(gen.dat, valid.dat, splr.dat))
      print("......")
      #--------------------------------------------------------------------------------
      results = rbind(results, compare_softImpute_splr_Kfold(gen.dat, valid.dat, splr.dat))
      print(".......")
      #--------------------------------------------------------------------------------
      results$true_rank = gen.dat$rank
      results$dim = dim[i]
      results$k = ncovariates
      #---------------------------------------------------------------------------
      # saving plots to disk
      #final.results = rbind(final.results, results)
      print(results)
      return(results)
   }
   if(ncores > 1)
      stopImplicitCluster()
   
   final.results$missinginess = missingness
   final.results$collinearity = coll
   
   if(missingness == 0){
      filename = "Compare_MC_Models_Mao_Simulation.csv"
   }else
      filename = paste0("Compare_MC_Models_Youssef_Simulation_",round(missingness*100),
                        "_coll_",coll, ".csv")
   
   write.csv(final.results, file = paste0(data_dir, filename), row.names = FALSE)
   print(paste("Results saved to",paste0(data_dir, filename)))
   return(final.results)
   #----------------------------
}

alpha_grid = c(1)
ncores = 6
lambda.1_grid = seq(0,2,length=20)
lambda.2_grid = seq(.9, 0, length=20) 


compare_and_save(0.9, TRUE,
                 lambda.1_grid = lambda.1_grid,lambda.2_grid = lambda.2_grid,
                 alpha_grid = alpha_grid, ncores=ncores)


compare_and_save(0.9, FALSE, 
                 lambda.1_grid = lambda.1_grid, lambda.2_grid = lambda.2_grid,
                 alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0.8, TRUE, 
                 lambda.1_grid = lambda.1_grid, lambda.2_grid = lambda.2_grid,
                 alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0.8, FALSE, 
                  lambda.1_grid = lambda.1_grid,lambda.2_grid = lambda.2_grid,
                  alpha_grid = alpha_grid, ncores=ncores)
compare_and_save(0,  FALSE,
                 lambda.1_grid = lambda.1_grid, lambda.2_grid = lambda.2_grid, 
                 alpha_grid = alpha_grid, ncores=ncores)

#------------------------------------------------------------------------------------  