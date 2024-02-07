simpute.cov.cv_splr <- function(Y, svdH, Y.valid, W_valid, lambda.factor=1/4, lambda.init=NA, n.lambda=20,
                           trace=FALSE, print.best=TRUE, tol=5, thresh=1e-5,
                           rank.init=10, rank.limit=50, rank.step=2, patience=2,
                            warm=NULL, quiet=FALSE){
   
   
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov(Y, X) * lambda.factor, lambda.init) 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   #----------------------------------------------------
   Y[Y==0] = NA
   m = dim(Y)[2]
   W_valid[W_valid ==1] = NA
   #W_valid[W_valid==0] =  1
   ys <- as(Y, "Incomplete")
   W_valid <- as(W_valid, "Incomplete")
   irow = W_valid@i
   pcol = W_valid@p
   W_valid = NULL
   #-----------------------------------------------------------------------
   rank.max <- rank.init
   best_fit <- list(error=Inf, rank_A=NA, rank_B=NA, lambda=NA, rank.max=NA)
   counter <- 1
   #---------------------------------------------------------------------
   for(i in seq(along=lamseq)) {
      fiti <- simpute.als.fit_splr(y = ys, svdH = svdH,
                                   trace=F, J=rank.max, thresh=thresh, lambda=lamseq[i],
                                   warm.start = warm, patience=patience, maxit=100)
                                   
      # get test estimates and test error
      M_valid = suvC(as.matrix(fiti$u),as.matrix(UD(fiti$v,fiti$d,m)),irow,pcol)
      #M_valid = (fiti$u %*% (fiti$d * t(fiti$v)))[W_valid==0]
      err = test_error(M_valid, Y_valid)
      rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
      #----------------------------
      warm <- fiti # warm start for next 
      if(trace==TRUE)
         print(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error on M = %.5f\n",
                       i, lamseq[i], rank.max, rank, err))
      #-------------------------
      # register best fir
      if(err < best_fit$error){
         best_fit$error = err
         best_fit$rank_B = rank
         best_fit$lambda = lamseq[i]
         best_fit$rank.max = rank.max
         best_fit$B_hat = M
         best_fit$best_fit = fiti
         best_fit$iter = i
         counter=1
      }else counter = counter + 1
      if(counter >= tol){
         if(quiet == FALSE)
            print(sprintf("Performance didn't improve for the last %d iterations.", counter))
         break
      }
      # compute rank.max for next iteration
      rank.max <- min(rank+rank.step, rank.limit)
   }
   #if(print.best==TRUE) print(best_fit)
   return(best_fit)
}
#---
simpute.cov.cv_splr_no_patience <- function(Y, svdH, Y.valid, W_valid, lambda.factor=1/4, lambda.init=NA, n.lambda=20,
                                trace=FALSE, print.best=TRUE, tol=5, thresh=1e-5,
                                rank.init=10, rank.limit=50, rank.step=2,
                                warm=NULL, quiet=FALSE,patience=2){
   
   
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov(Y, X) * lambda.factor, lambda.init) 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   #----------------------------------------------------
   Y[Y==0] = NA
   m = dim(Y)[2]
   W_valid[W_valid ==1] = NA
   Y <- as(Y, "Incomplete")
   W_valid <- as(W_valid, "Incomplete")
   irow = W_valid@i
   pcol = W_valid@p
   W_valid = NULL
   #-----------------------------------------------------------------------
   rank.max <- rank.init
   best_fit <- list()
   old_error = Inf
   old_fit = NULL
   #---------------------------------------------------------------------
   for(i in seq(along=lamseq)) {
      new_fit <- simpute.als.fit_splr(y = U, svdH = svdH,
                                   trace=F, J=rank.max, thresh=thresh, lambda=lamseq[i],
                                   warm.start = old_fit, patience=patience, maxit=100)
      
      M_valid = suvC(as.matrix(new_fit$u),as.matrix(UD(new_fit$v,new_fit$d,m)),irow,pcol)
      err = test_error(M_valid, Y_valid)
      #----------------------------
      if(err > old_error){
         best_fit$error = old_error
         best_fit$rank_B = rank
         best_fit$rank.max = rank.max
         best_fit$lambda = lamseq[i-1]
         best_fit$best_fit = old_fit
         best_fit$iter = i-1
         break
      }
      #-----------------------------------------
      rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
      old_fit <- new_fit # warm start for next 
      old_error = err
      #-------------------------------------------------------------------
      if(trace==TRUE)
         print(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error on M = %.5f\n",
                       i, lamseq[i], rank.max, rank, err))
      #-------------------------
      rank.max <- min(rank+rank.step, rank.limit)
      #----------------------------------------------------------------
      
   }
   #if(print.best==TRUE) print(best_fit)
   return(best_fit)
}