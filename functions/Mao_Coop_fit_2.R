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
   #list(A_hat = sout$A_hat)
   sout$A_hat
}

coop_fit <- function(A, X, Z, W, maxiter=100, epsilon=1e-3, 
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
                     trace_fin = FALSE, rank.step=4){
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
   Y = A * W
   W_valid <- matrix.split.train.test(W, testp=0.2)
   Y_train <-Y * W_valid
   Y_valid <- Y[W_valid==0]
   
   Y_tr = t(Y)
   W_tr = t(W)
   Y_valid_tr = t(Y_valid)
   Y_train_tr = t(Y_train)
   W_valid_tr = t(W_valid)
   yobs = Y_train != 0
   y.star <- Y_train
   Y_train_obs <- Y_train[yobs]
   #A_tr = t(A)
   
   
   # to obtain initial estimates, we run the MAO function on each covariate matrix individually
   A.hat_x = coop_fit_step(Y_train, X, W_valid, Y_valid, lambda1.grid, trace, rank.limit,
                           print.best, trace_fin, rank.step, tol=tol)
   
   # best_param = Mao.cv(A, X, Y, W, n_folds, lambda.1_grid, lambda.2_grid, alpha_grid, seed,
   #                       numCores, n1n2_optimized, theta_estimator)$best_parameters
   # A.hat_x = Mao.fit(Y, X, W, best_param$lambda.1, best_param$lambda.2, best_param$alpha, 
   #                   n1n2_optimized, theta_estimator )$A_hat
   
   
   A.hat_z = t(coop_fit_step(Y_train_tr, Z, W_valid_tr, Y_valid_tr, lambda1.grid, trace, rank.limit,
                           print.best, trace_fin, rank.step, tol=tol))
   # best_param = Mao.cv(A_tr, Z, Y_tr, W_tr, n_folds, lambda.1_grid, lambda.2_grid, alpha_grid, seed,
   #                       numCores, n1n2_optimized, theta_estimator)$best_parameters
   # A.hat_z = t(Mao.fit(Y_tr, Z, W_tr, best_param$lambda.1, best_param$lambda.2, best_param$alpha, 
   #                     n1n2_optimized, theta_estimator )$A_hat)
   
   # reporting test errors
   valid_error_x = test_error(A.hat_x[W_valid==0], Y_valid)
   valid_error_z = test_error(A.hat_z[W_valid==0], Y_valid)
   #print(test_error_z)
   valid_error_avg = test_error(((A.hat_x+A.hat_z)/2)[W_valid==0], Y[W_valid==0])
   print(paste("Validation error using only X is:",
         round(valid_error_x,5),
         ", and using only Z:",
         round(valid_error_z,5),
         ", and using their average:",
         round(valid_error_avg,5)
   ))
   
   # We assume to row covariates are minimized first for now.
   iter = 0
   n1 = dim(Y)[1]
   n2 = dim(Y)[2]
   A.hat = A.hat_x + A.hat_z
   rho.1 = 1/(1+rho)
   rho.2 = (1-rho)/(1+rho)
   part1 = rho.1 * Y_train_obs
   diff = Inf
   tol_counter = 0
   old_diff = Inf
   y.hat = Y
   y.hat = y.hat.old = Y#= (A.hat_x+A.hat_z)/2
   #yfill[ymiss] = y.hat[ymiss]
   # If X is better fit then we start with it, if Z is better then we start with Z
   # The two while loops are identical except for the order of the two models.
   
   # case 1: X is better
   valid_error_z = Inf # remove later!!!
   if(valid_error_x <= valid_error_z){
      print("X is better, Entering loop 1")
      while(iter < maxiter & diff > epsilon){
         # 1. row covariates, (Y1)   
         #y.star[yobs] = part1 - rho.2 * A.hat_z[yobs]  #* W * W_valid
         y.star = rho.1 * Y_train - rho.2 * A.hat_z  * W * W_valid
         
         A.hat_x = coop_fit_step(y.star, X, W_valid, Y_valid, lambda1.grid, trace, rank.limit,
                                 print.best, trace_fin, rank.step, tol=tol)
         
         # 2. Column Covariates
         y.star = rho.1 * Y_train - rho.2 * A.hat_x  * W * W_valid
         A.hat_z = t(coop_fit_step(t(y.star), Z, W_valid_tr, Y_valid_tr, lambda1.grid, trace, rank.limit,
                                   print.best, trace_fin, rank.step, tol=tol))
         
         
         
         
         y.hat = A.hat_x + A.hat_z
         diff = sqrt(mean((y.hat - y.hat.old)**2))
         y.hat.old = y.hat
         iter = iter + 1
         valid_error = test_error(y.hat[W_valid==0], Y_valid)
         test_error = test_error(y.hat[W==0], A[W==0])
         
         print(sprintf("Iteration %.0f - test error: %.3f - validation error: %.3f - diff: %.5f",
                       iter, test_error, valid_error, diff))
         # print(paste("Iteration",iter, "- test error =",test_error,
         #             "validation error =",valid_error, "- diff =",diff))
         
         if(diff >= old_diff){ tol_counter = tol_counter + 1}else{tol_counter = 0}
         if(tol_counter >= tol) break
         old_diff = diff
         # diff = sqrt(mean((A.hat - A.hat_x - A.hat_z)**2))
         # A.hat = A.hat_x + A.hat_z
         # A.hat_alt = (A_estim_x + A_estim_z) / 2
         # test_error_alt = test_error(A.hat_alt[W==0], A[W==0])
         
         # STOP if the difference has been increasing for 3 continuous iterations
      }
   }else { # case 2: Z is better
      print("Z is better, Entering loop 2")
      while(iter < maxiter & diff > epsilon){
         
         # 1. Column covariates (Y2)
         y.star = (rho.1 * Y_tr - rho.2 * t(A.hat_x)) 
         best_param = Mao.cv(A_tr, Z, y.star, W_tr, n_folds, lambda.1_grid, lambda.2_grid, alpha_grid, seed,
                             numCores, n1n2_optimized, theta_estimator)$best_parameters
         A.hat_z = t(Mao.fit(y.star* W_tr, Z, W_tr, best_param$lambda.1, best_param$lambda.2, best_param$alpha, 
                             n1n2_optimized, theta_estimator )$A_hat)
         # the following estimates A using Z according to the method in section "Estimation"
         A_estim_z = (1-rho) * A.hat_z + (1+rho) * A.hat_x
         # 2. row covariates, (Y1)   
         y.star = (rho.1 * Y - rho.2 * A.hat_z) 
         best_param = Mao.cv(A, X, y.star, W, n_folds, lambda.1_grid, lambda.2_grid, alpha_grid, seed,
                             numCores, n1n2_optimized, theta_estimator)$best_parameters
         A.hat_x = Mao.fit(y.star* W, X, W, best_param$lambda.1, best_param$lambda.2, best_param$alpha, 
                           n1n2_optimized, theta_estimator )$A_hat
         # the following estimates A using X according to the method in section "Estimation"
         A_estim_x = (1-rho) * A.hat_x + (1+rho) * A.hat_z
         # update stopping criteria
         diff = sqrt(mean((A.hat - A.hat_x - A.hat_z)**2))
         A.hat = A.hat_x + A.hat_z
         iter = iter + 1
         A.hat_alt = (A_estim_x + A_estim_z) / 2
         test_error_alt = test_error(A.hat_alt[W==0], A[W==0])
         test_error_Coop = test_error(A.hat[W==0], A[W==0])
         print(paste("Iteration",iter, "- test Error =",test_error_Coop,
                     "Error Alt =",test_error_alt, "- diff =",diff))
         
         # STOP if the difference has been increasing for 3 continuous iterations
         if(diff >= old_diff){ tol_counter = tol_counter + 1}else{tol_counter = 0}
         if(tol_counter >= tol) break
         old_diff = diff
         
      }
   }
   
   if(tol_counter>=tol) 
      print(paste0("Diff did not decrease for ", tol, " iterations. Quitting before reaching Max iterations."))
   
   if(iter >= maxiter){
      print("Reached Max iterations before converging.")
   }else if(tol_counter < tol)
      print("Converged.")
   
   test_error_Coop = test_error(A.hat[W==0], A[W==0])
   print(paste("Test error resulted from applying the Coop Model is",test_error_Coop))
   
   time_taken <- round(as.numeric(difftime(Sys.time(), start_time, units='mins')),3)   
   
   print(paste0("Run time is ",time_taken, " minutes."))
   print(paste0("Average Run time per iteration is ",round(time_taken/iter,3), " minutes."))
   return(A.hat)
} 


gen.dat <- generate_simulation_data_ysf(1,500,500,5,10, missing_prob = 0.8,coll=T, seed=2023)

out <- coop_fit(gen.dat$A, gen.dat$X, gen.dat$Z, gen.dat$W, rho=0.9, tol=2,trace_fin = F)






