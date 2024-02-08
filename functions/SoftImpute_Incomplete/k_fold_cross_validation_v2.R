simpute.cov.Kf_splr_no_patience_v2 <- function(Y, svdH, W, n_folds=5, lambda.factor=1/4, lambda.init=NA, n.lambda=20,
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
      W_fold[W == 0] = 1
      Y_train = Y * W_fold 
      #ymiss = W_fold==0
      
      #Y_train[Y_train==0] = NA
      # train_mask = Y_train
      # train_mask[train_mask != 0] = 1
      # train_mask = as(train_mask, "Incomplete")
      #Y_train = as(Y_train, "Incomplete")
      
      #W_fold[!ymiss] = NA
      #W_fold = as(W_fold, "Incomplete")
      #print(sum(W_fold==0))
      #print(length(Y[W_fold==0]))
      list(Y_train = Y_train, Y_valid = Y[W_fold==0], #irow=W_fold@i, pcol=W_fold@p,
           W_fold=W_fold)#, xbeta.obs=train_mask)# train_mask=train_mask, ymiss=ymiss)
   })
   #---------------------------------------------------------------------------
   # m = dim(Y)[2]
   # ymiss = Y==0
    Y[Y==0] = NA
    Y <- as(Y, "Incomplete")
   # W <- as(W, "Incomplete")
   # xbeta.obs <- M_estim <- Y
   #---------------------------------------------------------------------------
   rank.max <- rank.init
   best_fit <- list()
   old_error = Inf
   fiti = warm
   old_rank = NA
   #--------------------------------------------------------------------
   best_fit = list(best_fit=NULL)
   lambda = rank = rank.max = 0
   #---------------------------------------------------------------------
   for(fold in 1:n_folds){
      d = fold_data[[fold]]
      
      best_fit = simpute.cov.cv_splr_no_patience(d$Y_train, svdH, d$Y_valid,
                                                 d$W_fold, 
                                                 warm = best_fit$best_fit,
                                                 rank.init = rank.init,
                                                 trace = trace, rank.limit=rank.limit,
                                                 rank.step=rank.step,patience = patience)
      
      # d$xbeta.obs@x = best_fit$best_fit$xbeta.obs
      # if(fold>1){
      #    xbeta.obs.new = merge_sparse_mat(as.matrix(d$xbeta.obs), as.matrix(xbeta.obs.new))
      # }else
      #    xbeta.obs.new = as.matrix(d$xbeta.obs)
      
      lambda = lambda + best_fit$lambda
      rank.max = rank.max + best_fit$rank.max
      rank = rank + best_fit$rank_B
      best_fit$best_fit$xbeta.obs = NA
   }
   best_fit$lambda = lambda / n_folds
   best_fit$rank_B = as.integer(rank / n_folds)
   best_fit$rank.max = as.integer(rank.max / n_folds)
   
   # xbeta.obs.new[ymiss] = NA
   # xbeta.obs.new <- as(xbeta.obs.new, "Incomplete")
   # best_fit$best_fit$xbeta.obs <- xbeta.obs.new@x
   # best_fit$xbeta.obs = xbeta.obs.new
   # #print("Hi8")
   # #print(best_fit)
   # best_fit$best_fit <- simpute.als.fit_splr(y = Y, svdH = svdH,
   #                                    trace=trace, J=best_fit$rank_B, thresh=thresh, lambda=best_fit$lambda,
   #                                    warm.start = best_fit$best_fit, patience=patience, maxit=100)
   return(best_fit)
}

