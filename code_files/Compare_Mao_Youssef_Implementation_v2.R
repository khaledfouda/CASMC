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
                          trace=F, rank.limit = 30, quiet = TRUE,
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

compare_softImpute_L2 <- function(gen.dat, valid.dat, lambda.1_grid, rank.step, rank.limit,
                                  n.lambda){
   
   start_time = Sys.time()
   sout <- simpute.cov.cv(valid.dat$Y_train, gen.dat$X, valid.dat$W_valid,
                          valid.dat$Y_valid, trace=FALSE, rank.limit = rank.limit, 
                          print.best=FALSE, rank.step=rank.step, type="als", lambda1=0,
                          tol=2, n.lambda = n.lambda, quiet = TRUE)
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

compare_softImpute_Kfold <- function(gen.dat, lambda.1_grid, n_folds, rank.step, rank.limit,
                                     n.lambda){
   
   start_time = Sys.time()
   sout <- simpute.cov.kfold(gen.dat$Y, gen.dat$X, gen.dat$W, n_folds = n_folds, print.best = FALSE,
                             trace=FALSE, rank.limit = rank.limit, lambda1=0,n1n2 = 1,
                             warm=NULL,tol = 2, rank.step = rank.step, n.lambda = n.lambda)
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




compare_softImpute_splr <- function(gen.dat, valid.dat, splr.dat, rank.step, rank.limit, n.lambda){
   
   start_time = Sys.time()
   
   
   best_fit = simpute.cov.cv_splr(valid.dat$Y_train, splr.dat, valid.dat$Y_valid,
                                  valid.dat$W_valid, gen.dat$Y, trace=F, thresh=1e-6,
                                  n.lambda = n.lambda, rank.limit=rank.limit, maxit = 200,
                                  rank.step = rank.step, print.best = FALSE)
   
   
   fit1 = best_fit$fit1
   fit2 = best_fit$fit2
   sout = best_fit
   # get estimates and validate
   sout$B_hat = fit1$u %*% (fit1$d * t(fit1$v))
   sout$beta_hat =  fit2$u %*% (fit2$d * t(fit2$v))
   sout$A_hat = sout$B_hat + splr.dat$X %*% t(sout$beta_hat)
   
   results = list(model = "SoftImpute_splr")
   results$time =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
   results$alpha = NA
   results$lambda.1 = NA
   results$lambda.2 = sout$lambda
   results$error.test = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$error.all = test_error(sout$A_hat, gen.dat$A)
   results$error.beta = test_error(t(sout$beta_hat), gen.dat$beta)
   results$error.B = test_error(sout$B_hat, gen.dat$B)
   results$rank = qr(sout$A_hat)$rank
   results
}

compare_softImpute_splr_Kfold <- function(gen.dat, splr.dat, n_folds, 
                                          rank.step, rank.limit, n.lambda){
   
   start_time = Sys.time()
   
   best_fit = simpute.cov.Kf_splr(gen.dat$Y, splr.dat, gen.dat$W, trace = F,
                                  print.best = F, rank.limit=rank.limit,
                                  n.lambda = n.lambda, maxit = 200,
                                  n_folds = n_folds, rank.step=rank.step)
   
   
   fit1 = best_fit$fit1
   fit2 = best_fit$fit2
   sout = best_fit
   # get estimates and validate
   sout$B_hat = fit1$u %*% (fit1$d * t(fit1$v))
   sout$beta_hat =  fit2$u %*% (fit2$d * t(fit2$v))
   sout$A_hat = sout$B_hat + splr.dat$X %*% t(sout$beta_hat)   
   
   
   results = list(model = "SoftImpute_splr_Kfold")
   results$time =round(as.numeric(difftime(Sys.time(), start_time,units = "secs")))
   results$alpha = NA
   results$lambda.1 = NA
   results$lambda.2 = sout$lambda
   results$error.test = test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   results$error.all = test_error(sout$A_hat, gen.dat$A)
   results$error.beta = test_error(t(sout$beta_hat), gen.dat$beta)
   results$error.B = test_error(sout$B_hat, gen.dat$B)
   results$rank = qr(sout$A_hat)$rank
   results
}



compare_and_save <- function(missingness,coll=TRUE,  n_folds=3,
                             dim = seq(400,1000,200), ncovariates=10,
                             lambda.1_grid = seq(0,3,length=20),
                             lambda.2_grid = seq(.9, 0, length=20),
                             alpha_grid = seq(0.992, 1, length=10),ncores=1,
                             n.lambda=30, rank.limit=20, rank.step=2,
                             error_function = RMSE, seed=2024){
   
   
   ncores = min(ncores, length(dim))
   
   if (ncores > 1) {
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
   } else
      registerDoSEQ()
   
   print(paste("Running on",ncores,"core(s)."))
   test_error <<- error_function
   data_dir = "./saved_data/"
   stopifnot(missingness %in% c(0,0.8, 0.9))
   final.results = data.frame()
   
   for(i in 1:length(dim)){
   
   #final.results <- foreach(i = 1:length(dim), .combine='rbind') %do%  {
      set.seed(seed)
      results <- data.frame()
      if(missingness == 0){
         gen.dat <- generate_simulation_data_mao(n1=dim[i],n2=dim[i],m=ncovariates,r=ncovariates, seed=seed)
      }else
         gen.dat <- generate_simulation_data_ysf(2,dim[i],dim[i],ncovariates,ncovariates, 
                                                 missing_prob = missingness,coll=coll,seed=seed)
      #-------------------------------------------------------------------------------------
      # validation set to be used for the next two models
      valid.dat = list()
      valid.dat$W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2,seed = seed)
      valid.dat$Y_train <- gen.dat$Y * valid.dat$W_valid
      valid.dat$Y_valid <- gen.dat$Y[valid.dat$W_valid==0]
      #---------------------------------------------------------
      # SPLR data
      splr.dat = reduced_hat_decomp(gen.dat$X, 1e-2)
      gen.dat$X <- splr.dat$X
      #----------------------------------------------------------------------
      # fit 1. Mao
      print(i)
      results = rbind(results, compare_Mao(gen.dat, lambda.1_grid, lambda.2_grid,
                                           alpha_grid,1,n_folds))
      cat(" - M1 - ")
      #----------------------------------------------------------
      # soft Impute model without covariates
      results = rbind(results, compare_softImpute_orig(gen.dat, valid.dat))
      cat(" - M2 - ")
      #----------------------------------------------------------------------------
      # soft Impute model with covariates
      results = rbind(results, compare_softImpute_cov(gen.dat, valid.dat)) 
      cat(" - M3 - ")
      #-------------------------------------------------------------------------------------
      # Soft Impute with Covariates and With L2 regularization on the covariates
      results = rbind(results, compare_softImpute_L2(gen.dat, valid.dat, lambda.1_grid,
                                                     rank.step,rank.limit,n.lambda))
      cat(" - M4 -")
      #-------------------------------------------------------------------------------------
      # Soft Impute with Covariates and With L2 regularization on the covariates and K-fold cross-validation
      results = rbind(results, compare_softImpute_Kfold(gen.dat, lambda.1_grid, n_folds,
                                                       rank.step, rank.limit,
                                                       n.lambda))
      cat(" - M5 - ")
      #--------------------------------------------------------------------------------
      results = rbind(results, compare_softImpute_splr(gen.dat, valid.dat, splr.dat,
                                                       rank.step,rank.limit,n.lambda))
      cat(" - M6 - ")
      #--------------------------------------------------------------------------------
      results = rbind(results, compare_softImpute_splr_Kfold(gen.dat, splr.dat,
                                                             n_folds,rank.step, rank.limit,
                                                             n.lambda))
      cat(" - M7 -\n")
      #--------------------------------------------------------------------------------
      results$true_rank = gen.dat$rank
      results$dim = dim[i]
      results$k = ncovariates
      #---------------------------------------------------------------------------
      # saving plots to disk
      print(results)
      #return(results)
      final.results = rbind(final.results, results)
   }
   print("Exiting Loop ...")
   if(ncores > 1)
      stopCluster(cl)
      #stopImplicitCluster()
   
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

setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R", local=FALSE)

alpha_grid = c(1)
ncores = 1
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