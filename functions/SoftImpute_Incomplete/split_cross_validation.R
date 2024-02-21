simpute.cov.cv_splr <- function(Y, X_r, Y_valid, W_valid, lambda.factor=1/4, lambda.init=NA, n.lambda=20,
                           trace=FALSE, print.best=TRUE, tol=1, thresh=1e-6,
                           rank.init=10, rank.limit=50, rank.step=2,
                            warm=NULL, quiet=FALSE){
   
   
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov_splr(Y, X_r$svdH) * lambda.factor, lambda.init) 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   #----------------------------------------------------
   yfill =  Y
   ynas = Y == 0
   Y[ynas] = NA
   m = dim(Y)[2]
   Y <- as(Y, "Incomplete")
   xbeta.sparse = Y
   
   valid_ind = W_valid == 0
   #W_valid[W_valid == 1] = NA
   #W_valid[W_valid == 0] =  1
   #W_valid <- as(W_valid, "Incomplete")
   #irow = W_valid@i
   #pcol = W_valid@p
   W_valid = NULL
   #-----------------------------------------------------------------------
   rank.max <- rank.init
   best_fit <- list(error=Inf, rank_A=NA, rank_B=NA, lambda=NA, rank.max=NA)
   counter <- 0
   time_it <- rep(0,4)
   #---------------------------------------------------------------------
   for(i in seq(along=lamseq)) {
      start_time <- Sys.time()
      fiti <-  simpute.als.fit_splr(y=Y, svdH=X_r$svdH,  trace=F, J=rank.max,
                                    thresh=thresh, lambda=lamseq[i], return_obj = F, init = "naive",
                                    final.svd = T,maxit = 300, warm.start = warm)
      time_it[1] = time_it[1] + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
      start_time <- Sys.time()
      M = fiti$u %*% (fiti$d * t(fiti$v))
      xbeta.sparse@x <- fiti$xbeta.obs
      yfill[ynas] <- M[ynas]
      warm_xbeta = X_r$svdH$u %*% (X_r$svdH$v %*% yfill) 
      if(i>1) warm_xbeta = (warm_xbeta +  xbeta.estim) / 2
      #warm_xbeta = naive_MC(as.matrix(xbeta.sparse))
      warm_xbeta = propack.svd(warm_xbeta, X_r$rank)
      #if(i>1) warm_xbeta = fitx
      
      time_it[2] = time_it[2] + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
      start_time <- Sys.time()
      fitx <- simpute.als.splr.fit.nocov.fixedJ(xbeta.sparse, X_r$rank, maxit=300, final.trim = F,
                                                warm.start = warm_xbeta, trace.it=F, return_obj = F)
      time_it[3] = time_it[3] + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
      start_time <- Sys.time()
      xbeta.estim = fitx$u %*% (fitx$d * t(fitx$v))
      #xbeta.valid = suvC(as.matrix(fitx$u),as.matrix(UD(fitx$v,fitx$d,m)),irow,pcol)
      A_valid = M[valid_ind] + xbeta.estim[valid_ind]
      err = test_error(A_valid, Y_valid)
      time_it[4] = time_it[4] + as.numeric(difftime(Sys.time(), start_time,units = "secs"))                            
      rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
      #----------------------------
      warm <- fiti # warm start for next 
      if(trace==TRUE)
         print(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d,%d,%d, valid error = %.5f\n",
                       i, lamseq[i], rank.max, rank, fitx$iter, fiti$n_iter, err))
      #-------------------------
      # register best fir
      if(err < best_fit$error){
         best_fit$error = err
         best_fit$rank_B = rank
         best_fit$lambda = lamseq[i]
         best_fit$rank.max = rank.max
         best_fit$B_hat = M
         #best_fit$best_fit = fiti
         best_fit$xbeta = fitx
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
   print(time_it)
   #if(print.best==TRUE) print(best_fit)
   return(best_fit)
}
#---
