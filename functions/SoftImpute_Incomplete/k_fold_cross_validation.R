simpute.cov.Kf_splr <- function(Y, X_r, W,Px, n_folds=5, lambda.factor=1/4, lambda.init=NA, n.lambda=20,
                                            trace=FALSE, thresh=1e-5, maxit=100,
                                            rank.init=10, rank.limit=50, rank.step=2,
                                            warm=NULL, tol=1, print.best=TRUE){
   
   
   #----------------------------------------------------
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov_splr(Y, X_r$svdH) * lambda.factor, lambda.init) 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   #-----------------------------------------------------------------------
   # prepare the folds
   folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, W)
   fold_data <- lapply(1:n_folds, function(i) {
      W_fold = folds[[i]]
      valid_ind = W_fold==0 & W==1
      Y_train = Y * W_fold
      yfill <- Y_train
      ynas = Y_train == 0
      Y_train[ynas] = NA
      Y_train = as(Y_train, "Incomplete")
      xbeta.sparse = Y_train
      Y_valid = Y[valid_ind]

      list(Y_train = Y_train, Y_valid = Y_valid, yfill=yfill, ynas=ynas, valid_ind=valid_ind, xbeta.sparse=xbeta.sparse)
   })
   #---------------------------------------------------------------------------
   yfill =  Y
   ynas = Y == 0
   Y[ynas] = NA
   m = dim(Y)[2]
   Y <- as(Y, "Incomplete")
   xbeta.sparse = Y
   #---------------------------------------------------------------------------
   rank.max <- rank.init
   xbeta.estim = NULL
   best_error <- Inf
   best_rank <- best_lambda <- NA
   counter = 0
   #---------------------------------------------------------------------
   for(i in seq(along=lamseq)) {
      # initial fit
      
      # need Y_train (sparse); xbeta.sparse; yfill; ynas
      fiti <-  simpute.als.fit_splr(y=Y, svdH=X_r$svdH,  trace=F, J=rank.max,
                                    thresh=thresh, lambda=lamseq[i], return_obj = F, init = "naive",
                                    final.svd = T,maxit = maxit, warm.start = warm, Px=Px)
      M = fiti$u %*% (fiti$d * t(fiti$v))
      xbeta.sparse@x <- fiti$xbeta.obs
      # yfill[ynas] <- M[ynas]
      warm_xbeta =  as.matrix(X_r$X %*% fiti$beta.obs)
      # warm_xbeta = X_r$svdH$u %*% (X_r$svdH$v %*% yfill) 
      if(! is.null(xbeta.estim)) warm_xbeta = (warm_xbeta +  xbeta.estim) / 2
      warm_xbeta = propack.svd(warm_xbeta, X_r$rank)
      
      fitx <- simpute.als.splr.fit.nocov.fixedJ(xbeta.sparse, X_r$rank, maxit=maxit, final.trim = F,
                                                warm.start = warm_xbeta, trace.it=F, return_obj = F)
      xbeta.estim = fitx$u %*% (fitx$d * t(fitx$v))
      #-----------------------------------------------
      
      err <- rank <- 0
      for(fold in 1:n_folds){
         data = fold_data[[fold]]
         #-----
         fiti <-  simpute.als.fit_splr(y=data$Y_train, svdH=X_r$svdH,  trace=F, J=rank.max,
                                       thresh=thresh, lambda=lamseq[i], return_obj = F, init = "naive",
                                       final.svd = T,maxit = maxit, warm.start = warm, Px=Px)
         M = fiti$u %*% (fiti$d * t(fiti$v))
         data$xbeta.sparse@x <- fiti$xbeta.obs
         #data$yfill[data$ynas] <- M[data$ynas]
         
         #warm_xbeta = X_r$svdH$u %*% (X_r$svdH$v %*% data$yfill)
         warm_xbeta =  as.matrix(X_r$X %*% fiti$beta.obs)
         warm_xbeta = (warm_xbeta +  xbeta.estim) / 2
         warm_xbeta = propack.svd(warm_xbeta, X_r$rank)
         
         fitx <- simpute.als.splr.fit.nocov.fixedJ(data$xbeta.sparse, X_r$rank, maxit=maxit, final.trim = F,
                                                   warm.start = warm_xbeta, trace.it=F, return_obj = F)
         xbeta.estim = fitx$u %*% (fitx$d * t(fitx$v))
         A_valid = M[data$valid_ind] + xbeta.estim[data$valid_ind]
         err = err + test_error(A_valid, data$Y_valid)
         rank <- rank + sum(round(fiti$d, 4) > 0) # number of positive sing.values
      }
      err = err / n_folds
      rank = as.integer(rank / n_folds)
      #------------------------------------------------
      if(trace==TRUE)
         print(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
                       i, lamseq[i], rank.max, rank, err))
      #-----------------------------------------------------------------------
      if(err < best_error){
         best_error = err
         best_lambda = lamseq[i]
         best_rank = rank.max
         counter = 0
         best_xbeta = xbeta.estim
      }else 
         counter = counter + 1
      if(counter >= tol){
         print(sprintf("Performance didn't improve for the last %d iterations.", counter))
         break
      }
      #-------------------------------------------------------------------
      rank.max <- min(rank+rank.step, rank.limit)
      #----------------------------------------------------------------
      
   }
   if(print.best==TRUE) print(sprintf("lambda=%9.5g, rank.max = %d, error = %.5f\n",
                                      best_lambda, best_rank, best_error))
   
   # one last fit!
   # need Y_train (sparse); xbeta.sparse; yfill; ynas
   fiti <-  simpute.als.fit_splr(y=Y, svdH=X_r$svdH,  trace=F, J=best_rank,
                                 thresh=thresh, lambda=best_lambda, return_obj = F, init = "naive",
                                 final.svd = T,maxit = maxit, warm.start = warm, Px=Px)
   M = fiti$u %*% (fiti$d * t(fiti$v))
   xbeta.sparse@x <- fiti$xbeta.obs
   # yfill[ynas] <- M[ynas]
   # warm_xbeta = X_r$svdH$u %*% (X_r$svdH$v %*% yfill) 
   warm_xbeta =  as.matrix(X_r$X %*% fiti$beta.obs)
   warm_xbeta = (warm_xbeta +  best_xbeta) / 2
   warm_xbeta = propack.svd(warm_xbeta, X_r$rank)
   
   fitx <- simpute.als.splr.fit.nocov.fixedJ(xbeta.sparse, X_r$rank, maxit=maxit, final.trim = F,
                                             warm.start = warm_xbeta, trace.it=F, return_obj = F)
   xbeta.estim = fitx$u %*% (fitx$d * t(fitx$v))
   
   results = list()
   results$lambda2 = best_lambda
   results$B_hat = M
   results$xbeta_hat = xbeta.estim
   results$A_hat = M + xbeta.estim
   results$rank_A = qr(results$A_hat)$rank
   results$J = best_rank
   return(results)
}

