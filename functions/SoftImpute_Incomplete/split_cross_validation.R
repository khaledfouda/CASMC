simpute.cov.cv_splr <- function(Y, X_r, Y_valid, W_valid, y=NULL, lambda.factor=1/4, lambda.init=NA, n.lambda=20,
                           trace=FALSE, print.best=TRUE, tol=1, thresh=1e-6, maxit=100,
                           rank.init=3, rank.limit=50, rank.step=2,
                            warm=NULL, quiet=FALSE){
   
   #' Input: Y training partially observed with missing set to 0
   #' X_r list
   #' Y_valid  vector of validation values
   #' W_valid matrix where W_valid = 0 for validation and 1 for train/test
   #' y (optional) whole training set including the validation for one last fit at the end.
   #' --------------------------------------------------------------------
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov_splr(Y, X_r$svdH) * lambda.factor, lambda.init) 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   #----------------------------------------------------
   Y[Y == 0] = NA
   #m = dim(Y)[2]
   Y <- as(Y, "Incomplete")
   xbeta.sparse = Y
   irow = Y@i
   pcol = Y@p
   
   #valid_ind = W_valid == 0
   W_valid[W_valid == 1] = NA
   W_valid[W_valid == 0] =  1
   W_valid <- as(W_valid, "Incomplete")
   virow = W_valid@i
   vpcol = W_valid@p
   W_valid = NULL
   #-------------
   # Initialize warm.start for the second model
   svdX = fast.svd(X_r$X)
   Ux = svdX$u
   Vx = svdX$d * t(svdX$v)
   X0 = ginv(t(Vx)%*%Vx) %*% t(Vx)
   warm.start.beta = list()
   warm.start.beta$X1 = X0 %*% t(Ux)
   warm.start.beta$X2 = X0 %*% Vx
   Xinv = ginv(X_r$X)
   
   #-----------------------------------------------------------------------
   rank.max <- rank.init
   best_fit <- list(error=Inf, rank_A=NA, rank_B=NA, lambda=NA, rank.max=NA)
   counter <- 0
   #---------------------------------------------------------------------
   for(i in seq(along=lamseq)) {
      fiti <-  simpute.als.fit_splr(y=Y, svdH=X_r$svdH,  trace=F, J=rank.max,
                                    thresh=thresh, lambda=lamseq[i], init = "naive",
                                    final.svd = T,maxit = maxit, warm.start = warm)
      xbeta.sparse@x <- fiti$xbeta.obs
      #---------
      # prepare warm.start.beta:
      if(i == 1){
         B = t( Xinv %*% naive_MC(as.matrix(xbeta.sparse))) # B = (X^-1 Y)'
         warm.start.beta$Bsvd = fast.svd(B)
      }else warm.start.beta$Bsvd = fitx
      #---------------------------
      # fit second model:
      fitx = simpute.als.splr.fit.beta(xbeta.sparse, X_r$X, X_r$rank, final.trim = F, thresh=thresh,
                                       warm.start = warm.start.beta, trace.it = F,maxit=maxit)
      #--------------------------------------------------------------
      # predicting validation set and xbetas for next fit:
      Xv = X_r$X %*% fitx$v
      fiti$xbeta.obs <- suvC(Xv, t(fitx$d * t(fitx$u)), irow, pcol)
      Xbvalid = suvC(Xv, t(fitx$d * t(fitx$u)), virow, vpcol)
      Mvalid = suvC(fiti$u, t(fiti$d * t(fiti$v)), virow, vpcol)
      #--------------------------------------------
      err = test_error(Mvalid+Xbvalid, Y_valid)
      rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
      
      
      #---------------------------------------------------------------------
      #----------------------------
      warm <- fiti # warm start for next 
      if(trace==TRUE)
         print(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d,%d, valid error = %.5f",
                       i, lamseq[i], rank.max, rank, fiti$n_iter, err))
      #-------------------------
      # register best fir
      if(err < best_fit$error){
         best_fit$error = err
         best_fit$rank_B = rank
         best_fit$lambda = lamseq[i]
         best_fit$rank.max = rank.max
         best_fit$fit1 = fiti
         best_fit$fit2 = fitx
         best_fit$iter = i
         counter=0
      }else counter = counter + 1
      if(counter >= tol){
         if(quiet == FALSE)
            print(sprintf("Performance didn't improve for the last %d iterations.", counter))
         break
      }
      # compute rank.max for next iteration
      rank.max <- min(rank+rank.step, rank.limit)
   }
   # fit one last time full model, if the train/valid is provided
   if(! is.null(y)){
      best_fit$fit1$xbeta.obs <- suvC(X_r$X %*% best_fit$fit2$v, 
                                      t(best_fit$fit2$d * t(best_fit$fit2$u)),
                                      y@i, y@p)
      best_fit$fit1 <-  simpute.als.fit_splr(y=y, svdH=X_r$svdH,  trace=F, J=best_fit$rank.max,
                                    thresh=thresh, lambda=best_fit$lambda, init = "naive",
                                    final.svd = T,maxit = maxit, warm.start = best_fit$fit1)
      xbeta.sparse = y
      xbeta.sparse@x <- best_fit$fit1$xbeta.obs
      #B = t( Xinv %*% naive_MC(as.matrix(xbeta.sparse))) # B = (X^-1 Y)'
      #warm.start.beta$Bsvd = fast.svd(B)
      warm.start.beta$Bsvd = fitx
      best_fit$fit2 = simpute.als.splr.fit.beta(xbeta.sparse, X_r$X, X_r$rank, final.trim = F, thresh=thresh,
                                       warm.start = warm.start.beta, trace.it = F,maxit=maxit)

   }
   return(best_fit)
}
#---
