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
   W <- W * W_valid
   Y <- Y * W 
   
   #Y_t <- t(Y)
   Y_valid_tr = t(Y_valid)
   W_valid_tr = t(W_valid)
   ymiss <- W_valid == 0
   
   # to obtain initial estimates, we run the MAO function on each covariate matrix individually
   A.hat_x = coop_fit_step(Y, X, W_valid, Y_valid, lambda1.grid, trace, rank.limit,
                           print.best, trace_fin, rank.step, tol=tol)

   A.hat_z = t(coop_fit_step(t(Y), Z, W_valid_tr, Y_valid_tr, lambda1.grid, trace, rank.limit,
                           print.best, trace_fin, rank.step, tol=tol))
   # reporting test errors
   valid_error_x = test_error(A.hat_x[ymiss], Y_valid)
   valid_error_z = test_error(A.hat_z[ymiss], Y_valid)
   #print(test_error_z)
   valid_error_avg = test_error(((A.hat_x+A.hat_z)/2)[ymiss], Y_valid)
   print(paste("Validation error using only X is:",
         round(valid_error_x,4),
         ", and using only Z:",
         round(valid_error_z,4),
         ", and using their average:",
         round(valid_error_avg,4)
   ))
   
   # We assume to row covariates are minimized first for now.
   #n1 = dim(Y)[1]
   #n2 = dim(Y)[2]
   #A.hat = A.hat_x + A.hat_z
   
   iter = 0
   rho.1 = 1/(1+rho)
   rho.2 = (1-rho)/(1+rho)
   Y = rho.1 * Y
   
   diff <- old_diff <- Inf
   tol_counter = 0
   y.hat.old = A.hat_x + A.hat_z
   X_first <- ifelse(valid_error_x <= valid_error_z, TRUE, FALSE)
   # If X is better fit then we start with it, if Z is better then we start with Z
   # The two while loops are identical except for the order of the two models.
   
   # case 1: X is better
   if(valid_error_x <= valid_error_z){
      print("X is better, Entering loop 1")
   }else
      print("Z is better, Entering loop 2")
   
   
   while(iter < maxiter & diff > epsilon){
      if(X_first){
         # 1. row covariates, (Y1)   
         y.star = Y - rho.2 * (A.hat_z * W) 
         A.hat_x = coop_fit_step(y.star, X, W_valid, Y_valid, lambda1.grid, trace, rank.limit,
                                 print.best, trace_fin, rank.step, tol=tol)
         # 2. Column Covariates
         y.star = Y - rho.2 * (A.hat_x * W) 
         A.hat_z = t(coop_fit_step(t(y.star), Z, W_valid_tr, Y_valid_tr, lambda1.grid, trace, rank.limit,
                                   print.best, trace_fin, rank.step, tol=tol))
      }else{
         # 1. Column Covariates (Z)
         y.star = Y - rho.2 * (A.hat_x * W) 
         A.hat_z = t(coop_fit_step(t(y.star), Z, W_valid_tr, Y_valid_tr, lambda1.grid, trace, rank.limit,
                                   print.best, trace_fin, rank.step, tol=tol))
         # 2. row covariates, (X)   
         y.star = Y - rho.2 * (A.hat_z * W) 
         A.hat_x = coop_fit_step(y.star, X, W_valid, Y_valid, lambda1.grid, trace, rank.limit,
                                 print.best, trace_fin, rank.step, tol=tol)
      }
      #-----------------------------------------------------
      # prediction is the sum of two
      y.hat = A.hat_x + A.hat_z
      # compute the difference in prediction
      diff = sqrt(mean((y.hat[ymiss] - y.hat.old[ymiss])**2))
      # implement the objective function later!!!
      #--------------------------------------------------------
      # update variables and print
      y.hat.old = y.hat
      iter = iter + 1
      valid_error = test_error(y.hat[ymiss], Y_valid)
      print(sprintf("Iteration %.0f - validation error: %.3f - diff: %.5f",
                    iter, valid_error, diff))
      if(diff >= old_diff){ tol_counter = tol_counter + 1}else{tol_counter = 0}
      if(tol_counter >= tol) break
      old_diff = diff
      # STOP if the difference has been increasing for 3 continuous iterations
   }
   
   if(tol_counter>=tol) 
      print(paste0("Diff did not decrease for ", tol, " iterations. Quitting before reaching Max iterations."))
   
   if(iter >= maxiter){
      print("Reached Max iterations before converging.")
   }else if(tol_counter < tol)
      print("Converged.")
   
   valid_error = test_error(y.hat[ymiss], Y_valid)
   print(paste("Validation error resulted from applying the Coop Model is",valid_error))
   
   time_taken <- round(as.numeric(difftime(Sys.time(), start_time, units='mins')),3)   
   
   print(paste0("Run time is ",time_taken, " minutes."))
   print(paste0("Average Run time per iteration is ",round(time_taken/iter,3), " minutes."))
   return(y.hat)
} 


gen.dat <- generate_simulation_data_ysf(1,500,500,5,5, missing_prob = 0.8,coll=F)
W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_valid <- gen.dat$Y[W_valid==0]

out <- coop_fit(gen.dat$Y, gen.dat$X, gen.dat$Z, gen.dat$W, W_valid, Y_valid, rho=0.9, tol=2,trace_fin = F)

test_error(out[gen.dat$W==0], gen.dat$A[gen.dat$W==0])

