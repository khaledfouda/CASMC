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
      train_mask = Y_train
      train_mask[!is.na(train_mask)] = 1
      train_mask = as(train_mask, "Incomplete")
      Y_train = as(Y_train, "Incomplete")
      
      W_fold[!ymiss] = NA
      W_fold = as(W_fold, "Incomplete")

      list(Y_train = Y_train, Y_valid = Y[ymiss], irow=W_fold@i, pcol=W_fold@p,
           W_fold=W_fold, xbeta.obs=train_mask, train_mask=train_mask, ymiss=ymiss)
   })
   #---------------------------------------------------------------------------
   m = dim(Y)[2]
   ymiss = Y==0
   Y[ymiss] = NA
   Y <- as(Y, "Incomplete")
   W <- as(W, "Incomplete")
   xbeta.obs <- M_estim <- Y
   #---------------------------------------------------------------------------
   rank.max <- rank.init
   best_fit <- list()
   old_error = Inf
   fiti = warm
   old_rank = NA
   #---------------------------------------------------------------------
   for(i in seq(along=lamseq)) {
      #print("Hi")
      M_estim@x = 100
      fiti <- simpute.als.fit_splr(y = Y, svdH = svdH,
                           trace=F, J=rank.max, thresh=thresh, lambda=lamseq[i],
                           warm.start = fiti, patience=patience, maxit=100)
      xbeta.obs@x = fiti$xbeta.obs
      #print("Hi2")
      err <- rank <- 0
      for(fold in 1:n_folds){
         data = fold_data[[fold]]
         #   fiti$xbeta.obs = data$xbeta.obs
         #if(i > 1){
            fiti$xbeta.obs = (xbeta.obs * data$train_mask)@x
         #}else
         #   fiti = NULL
         #print(range(fiti$xbeta.obs))
         #fiti$xbeta.obs = rep(0, length(fiti$xbeta.obs))
         # print(length(fiti$xbeta.obs))
         # print(length(data$train_mask@x))
         # print(fiti$xbeta.obs[1:10])
         # print(length(data$Y_train@x))
         #print("Hi3")
         #return(list(xbeta.obs, data$train_mask))
         fiti <- simpute.als.fit_splr(y = data$Y_train, svdH = svdH,
                                         trace=F, J=rank.max, thresh=thresh, lambda=lamseq[i],
                                         warm.start = fiti, patience=patience, maxit=300)
         #print("Hi4")
         data$xbeta.obs@x = fiti$xbeta.obs
         #print("Hi5")
         if(fold>1){
            #return(list(a=as.matrix(data$xbeta.obs), b= as.matrix(xbeta.obs.new)))
            xbeta.obs.new = merge_sparse_mat(as.matrix(data$xbeta.obs), as.matrix(xbeta.obs.new))
         }else
            xbeta.obs.new = as.matrix(data$xbeta.obs)
         #print("Hi6")
         
         M_valid = suvC(as.matrix(fiti$u),as.matrix(UD(fiti$v,fiti$d,m)),data$irow,data$pcol)
         M_estim[data$ymiss] = M_valid
         #err = err + test_error(M_valid, gen.dat$B[data$ymiss])#data$Y_valid)
         rank <- rank + sum(round(fiti$d, 4) > 0) # number of positive sing.values
      }
      #print("Hi7")
      err = test_error(M_estim@x , Y[!ymiss]) #data$Y_valid)
      #err = err / n_folds
      rank = as.integer(rank / n_folds)
      xbeta.obs.new[ymiss] = NA
      #return(xbeta.obs.new)
      #print(length(as(xbeta.obs.new, "Incomplete")@x))
      xbeta.obs.new <- as(xbeta.obs.new, "Incomplete")
       fiti$xbeta.obs <- xbeta.obs.new@x
      #------------------------------------------------
      if(trace==TRUE)
         print(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error on M = %.5f\n",
                       i, lamseq[i], rank.max, rank, err))
      #-----------------------------------------------------------------------
      if(err >= old_error| i==length(lamseq)){
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
   #print("Hi8")
   #print(best_fit)
   best_fit$best_fit <- simpute.als.fit_splr(y = Y, svdH = svdH,
                                      trace=T, J=best_fit$rank.max, thresh=thresh, lambda=best_fit$lambda,
                                      warm.start = fiti, patience=patience, maxit=100)
   return(best_fit)
}


merge_sparse_mat <- function(matrix1, matrix2){
   result <- matrix1
   #non_zero_both <- (!is.na(matrix1)) & (! is.na(matrix2))
   #result[non_zero_both] <- (matrix1[non_zero_both] + matrix2[non_zero_both]) / 2
   result[is.na(matrix1) & (!is.na(matrix2))] <- matrix2[is.na(matrix1) & (!is.na(matrix2))]
   result[is.na(matrix1) & (is.na(matrix2))] <- 0
   result
}
# 
# matrix1 = best_fit[[1]]
# matrix2 = best_fit[[2]]
# 
# 
# sum(is.na(matrix2))
