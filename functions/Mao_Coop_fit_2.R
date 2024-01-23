# Make sure to run import_lib.R before running this file.


coop_fit_step <- function(Y_train, X, W_valid, Y_valid, lambda1.grid = seq(0,20,length.out=20),
                          trace=FALSE, rank.limit=30, print.best=FALSE, trace_fin = FALSE,
                          rank.step=4, type="als", tol=2){
   
   sout <- simpute.cov.cv(Y_train, X, W_valid, Y_valid, trace=trace, rank.limit = rank.limit, 
                          print.best=print.best, rank.step=rank.step, type=type, lambda1=0, tol=tol,
                          quiet=TRUE)
   
   sout <- simpute.cov.cv.lambda1(Y_train, X, W_valid, Y_valid, sout$lambda, sout$rank.max, 
                                  print.best = print.best,
                                  trace=trace, lambda1.grid = lambda1.grid ,n1n2 = 1, warm=NULL)
   
   if(trace_fin)
      print(sprintf("Validation Error: %.3f, rank: %.0f",
                          test_error(sout$A_hat[W_valid==0], Y_valid), sout$rank_A))
   #list(A = sout$A_hat, beta = sout$beta_hat, B=sout$B_hat, lambda1 = sout$lambda1, lambda2 = sout$lambda2)
   #sout$A_hat
   sout
}

coop_fit <- function(Y, X, Z, W, W_valid, Y_valid, maxiter=100, epsilon=1e-3, 
                         rho = 1, # the agreement penalty
                         n_folds=5,
                         lambda1.grid = seq(0,2,length=10),
                         lambda2.grid = seq(.9, 0, length=10),
                         alpha_grid = c(1),   #seq(0.992, 1, length=10),
                         numCores = 4,
                         n1n2_optimized = TRUE,
                         theta_estimator = theta_default,
                         seed = 2023,
                         tol=3, trace=FALSE, rank.limit = 30, print.best=FALSE,
                     trace_fin = FALSE, rank.step=4,
                     early_stopping=TRUE, patience =3,
                     verbose=TRUE, track=TRUE){
   #' Input:
   #'       A: The response matrix of dimension n1 by n2. We assume that Y = A * W
   #'       X: The row covariates of dimension n1 by m1
   #'       Z: The column covariates of dimension n2 by m2
   #'       W: Is the mask matrix with the same dimension as Y
   #'       niter:
   #' 
   #' Output:
   #'       A.hat
   #' ----------------------------------------------------------------
   # The following is to compute the time taken to apply the algorithm till convergence or reaching max iter
   start_time <- Sys.time()
   # Compute some variables and have the transposes ready.
   A = Y
   W <- W * W_valid
   Y <- Y * W 
   
   #Y_t <- t(Y)
   Y_valid_tr = t(Y_valid)
   W_valid_tr = t(W_valid)
   ymiss <- W_valid == 0
   
   # to obtain initial estimates, we run the MAO function on each covariate matrix individually
   fit.x = coop_fit_step(Y, X, W_valid, Y_valid, lambda1.grid, trace, rank.limit,
                           print.best, trace_fin, rank.step, tol=tol)
   A.hat_x = fit.x$A_hat
   fit.z = coop_fit_step(t(Y), Z, W_valid_tr, Y_valid_tr, lambda1.grid, trace, rank.limit,
                           print.best, trace_fin, rank.step, tol=tol)
   A.hat_z = t(fit.z$A_hat)
   y.hat <- y.hat.old <- ((A.hat_x+A.hat_z)/2)
   # reporting test errors
   valid_error_x = test_error(A.hat_x[ymiss], Y_valid)
   valid_error_z = test_error(A.hat_z[ymiss], Y_valid)
   valid_error_avg = test_error(y.hat[ymiss], Y_valid)
   
   if(verbose)
      print(paste("Validation error using only X is:",
            round(valid_error_x,4),
            ", and using only Z:",
            round(valid_error_z,4),
            ", and using their average:",
            round(valid_error_avg,4)
      ))
   
   best_fit = list(fit_x = fit.x, fit_z = fit.z, preds=y.hat, valid_error=valid_error_avg)
   #-------------------------------
   iter = 1
   rho.1 = 1/(1+rho)
   rho.2 = (1-rho)/(1+rho)
   Y = rho.1 * Y
   
   diff <- old_diff <- Inf
   tol_counter = 0
   X_first <- ifelse(valid_error_x <= valid_error_z, TRUE, FALSE)
   # If X is better fit then we start with it, if Z is better then we start with Z
   # The two while loops are identical except for the order of the two models.
   
   
   while(iter < maxiter & diff > epsilon){
      if(!X_first){
         # 1. row covariates, (Y1)   
         y.star = Y - rho.2 * (A.hat_z * W) 
         fit.x = coop_fit_step(y.star, X, W_valid, Y_valid, lambda1.grid, trace, rank.limit,
                                 print.best, trace_fin, rank.step, tol=tol)
         A.hat_x = fit.x$A_hat
         # 2. Column Covariates
         y.star = Y - rho.2 * (A.hat_x * W) 
         fit.z = coop_fit_step(t(y.star), Z, W_valid_tr, Y_valid_tr, lambda1.grid, trace, rank.limit,
                                   print.best, trace_fin, rank.step, tol=tol)
         A.hat_z = t(fit.z$A_hat)
      }else{
         # 2. Column Covariates
         y.star = Y - rho.2 * (A.hat_x * W) 
         fit.z = coop_fit_step(t(y.star), Z, W_valid_tr, Y_valid_tr, lambda1.grid, trace, rank.limit,
                               print.best, trace_fin, rank.step, tol=tol)
         A.hat_z = t(fit.z$A_hat)
         # 2. Column Covariates (Z)
         y.star = Y - rho.2 * (A.hat_z * W) 
         fit.x = coop_fit_step(y.star, X, W_valid, Y_valid, lambda1.grid, trace, rank.limit,
                               print.best, trace_fin, rank.step, tol=tol)
         A.hat_x = fit.x$A_hat
      }
      #-----------------------------------------------------
      # prediction is the sum of two
      y.hat = A.hat_x + A.hat_z
      # compute the difference in prediction
      diff = sqrt(mean((y.hat[ymiss] - y.hat.old[ymiss])**2))
      valid_error = test_error(y.hat[ymiss], Y_valid)
      
      if(track){
         
      # loss as defined in (Y0)
      loss <- norm( (A - A.hat_x - A.hat_z) * W , type="F")^2 +
         rho * norm( (A.hat_x - A.hat_z) * W, type="F") ^2 +
         fit.x$lambda1 * norm(fit.x$beta_hat, type="F")^2 +
         fit.z$lambda1 * norm(fit.z$beta_hat, type="F") ^2 +
         fit.x$lambda2 * sum(svd(fit.x$B_hat)$d) +
         fit.z$lambda2 * sum(svd(fit.z$B_hat)$d)
      
      print(sprintf("Iteration %.0f - validation error: %.4f - diff: %.3f - loss: %.3f",
                    iter, valid_error, diff, loss))
      }
      #--------------------------------------------------------
      # update variables and print
      if(valid_error < best_fit$valid_error){
         tol_counter = 0
         best_fit = list(fit_x = fit.x, fit_z = fit.z, preds=y.hat, valid_error=valid_error)
      }else 
         tol_counter = tol_counter + 1
      if(early_stopping){
         if(tol_counter >= patience)
            break
      }
      
      
      # update variables for next iteration
      iter = iter + 1
      y.hat.old = y.hat
      old_diff = diff
   }
   if(verbose){
      if(tol_counter>=patience) 
         print(paste0("Validation error has not decreased for ",
                      patience, " iterations. Quitting before reaching Max iterations."))
      if(iter >= maxiter){
         print("Reached Max iterations before converging.")
      }else if(tol_counter < patience)
         print(paste0("Converged after ", iter, " iterations."))
      print(paste("Validation error resulted from best fit is",round(best_fit$valid_error,4)))
      time_taken <- as.numeric(difftime(Sys.time(), start_time, units='mins'))   
      print(paste0("Total run time is ",round(time_taken,1), " minutes."))
      print(paste0("Average run time per iteration is ",round(time_taken/iter,1), " minutes."))
   }
   return(best_fit)
} 



coop_find_rho <- function(gen.dat, W_valid, early_maxiter=3, final_maxiter=30,
                          rho.grid = seq(0.1,0.99,length.out=10), print_best=TRUE,
                          seed=2023, patience=5, verbose=FALSE, tol=2, track=FALSE,
                          max_cores=10){
   num_cores = length(rho.grid)
   if(max_cores < num_cores){
      num_cores = round(length(rho.grid)/2)
      if(max_cores < num_cores){
         num_cores = min(round(length(rho.grid)/3),max_cores)
      }
   }
   print(paste("Running on",num_cores,"cores."))
   require(foreach)
   require(doParallel)
   
   Y_valid <- gen.dat$Y[W_valid==0]
   start_time <- Sys.time()
   
   cl <- makeCluster(num_cores)
   registerDoParallel(cl)
   
   #for(rho in rho.grid){
   results <- foreach(rho = rho.grid, .combine='list', .multicombine = TRUE,
                      .export = c("coop_fit", "coop_fit_step")) %dopar%{   
      #if(print_best)
      #   print(paste0("Attempting rho = ",round(rho,4)))
    source("./code_files/import_lib.R")
      out <- coop_fit(gen.dat$Y, gen.dat$X, gen.dat$Z, gen.dat$W, W_valid,
                      Y_valid, rho=rho, tol=tol,trace_fin = F, verbose=verbose,
                      early_stopping = TRUE, track=track,
                      patience=patience, maxiter = early_maxiter, seed = seed)
      list(rho=rho, fit=out)
   }
   #stopImplicitCluster()
   stopCluster(cl)
   #    if (out$valid_error < lowest_error){
   #       lowest_error <- out$valid_error
   #       best_rho <- rho
   # }
   lowest_error = Inf
   print("done")
   for(res in results){
      #print(res)
      if(res$fit$valid_error < lowest_error){
         lowest_error <- res$fit$valid_error
         best_rho = res$rho
         best_fit = res$fit
      }
   }
   
   if(print_best)
      print(paste("Best rho is",round(best_rho,2),
                   "with validation error of", round(lowest_error,4) ))
   # 
   # out <- coop_fit(gen.dat$Y, gen.dat$X, gen.dat$Z, gen.dat$W, W_valid,
   #                 Y_valid, rho=best_rho, tol=tol,trace_fin = F, verbose=print_best,
   #                 early_stopping = TRUE, track=print_best,
   #                 patience=patience, maxiter = final_maxiter, seed = seed)
   # 
   
   time_in_minutes =  as.numeric(difftime(Sys.time(), start_time, units='mins'))
   print(paste("Total time taken: ", round(time_in_minutes,1),"minutes."))
   list(rho=best_rho, fit=best_fit, time_in_minutes=time_in_minutes)
}




setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

gen.dat <- generate_simulation_data_ysf(1,500,500,5,10, missing_prob = 0.9,coll=T,seed=3023)
W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)

out <- coop_find_rho(gen.dat, W_valid,  print_best = TRUE,early_maxiter = 50,max_cores = 10,
                     rho.grid = seq(0.1,0.99,length.out=10))
out$time_in_minutes
out$rho
 
test_error(out$fit$preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0])

Y_valid <- gen.dat$Y[W_valid==0]

out <- coop_fit(gen.dat$Y, gen.dat$X, gen.dat$Z, gen.dat$W, W_valid,
                Y_valid, rho=0.3, tol=2,trace_fin = F, verbose=FALSE, early_stopping = TRUE,
                patience=5, maxiter = 15, seed = 3023)

test_error(out$preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
test_error(out$preds[W_valid==0], gen.dat$A[W_valid==0])

