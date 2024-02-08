simpute.cov.Kf_splr_no_patience <- function(Y, svdH, W, n_folds=5, lambda.factor=1/4, lambda.init=NA, n.lambda=20,
                                            trace=FALSE, thresh=1e-5,
                                            rank.init=10, rank.limit=50, rank.step=2,
                                            warm=NULL, patience=2){
   
   
   #----------------------------------------------------
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov_splr(Y, svdH) * lambda.factor, lambda.init) 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   #-----------------------------------------------------------------------
   # prepare the folds
   folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, W)
   fold_data <- lapply(1:n_folds, function(i) {
      W_fold = folds[[i]]
      ymiss = W_fold==0 & W==1
      
      Y_train = Y * W_fold
      Y_train[ymiss] = NA
      Y_train[W==0] = NA
      Y_train = as(Y_train, "Incomplete")
      
      W_fold[!ymiss] = NA
      W_fold = as(W_fold, "Incomplete")

      list(Y_train = Y_train, Y_valid = Y[ymiss], irow=W_fold@i, pcol=W_fold@p, xbeta.obs=NULL)
   })
   #---------------------------------------------------------------------------
   m = dim(Y)[2]
   Y[Y==0] = NA
   Y <- as(Y, "Incomplete")
   #---------------------------------------------------------------------------
   rank.max <- rank.init
   best_fit <- list()
   old_error = Inf
   fiti = warm
   old_rank = NA
   #---------------------------------------------------------------------
   for(i in seq(along=lamseq)) {
      
      err <- rank <- 0
      for(fold in 1:n_folds){
         data = fold_data[[fold]]
         #if(i > 1)
         #   fiti$xbeta.obs = data$xbeta.obs
         fiti <- simpute.als.fit_splr(y = data$Y_train, svdH = svdH,
                                         trace=F, J=rank.max, thresh=thresh, lambda=lamseq[i],
                                         warm.start = data$xbeta.obs, patience=patience, maxit=300)
         fold_data[[fold]]$xbeta.obs = fiti#$xbeta.obs
         M_valid = suvC(as.matrix(fiti$u),as.matrix(UD(fiti$v,fiti$d,m)),data$irow,data$pcol)
         err = err + test_error(M_valid, data$Y_valid)
         rank <- rank + sum(round(fiti$d, 4) > 0) # number of positive sing.values
      }
      err = err / n_folds
      rank = as.integer(rank / n_folds)
      #------------------------------------------------
      if(trace==TRUE)
         print(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error on M = %.5f\n",
                       i, lamseq[i], rank.max, rank, err))
      #-----------------------------------------------------------------------
      if(err >= old_error){
         best_fit$error = old_error
         best_fit$rank_B = old_rank
         best_fit$rank.max = rank.max
         best_fit$best_fit = fiti
         best_fit$lambda = lamseq[i-1]
         best_fit$iter = i - 1
         break
      }
      #-----------------------------------------
      old_error = err
      old_rank = rank
      #-------------------------------------------------------------------
      rank.max <- min(rank+rank.step, rank.limit)
      #----------------------------------------------------------------
      
   }
   #best_fit$best_fit <- simpute.als.fit_splr(y = Y, svdH = svdH,
   #                                     trace=F, J=best_fit$rank.max, thresh=thresh, lambda=best_fit$lambda,
   #                                     warm.start = NULL, patience=patience, maxit=100)
   return(best_fit)
}
