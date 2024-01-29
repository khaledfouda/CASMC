simpute.cov.cv_splr <- function(Y, X, W, Y.valid, lambda.factor=1/4, lambda.init=NA, n.lambda=20,
                           trace=FALSE, print.best=TRUE, tol=5, thresh=1e-5,
                           rank.init=10, rank.limit=50, rank.step=2, patience=2,
                           lambda1=0, n1n2=1, warm=NULL, quiet=FALSE){
   
   
   stopifnot(n1n2 %in% 1:3)
   # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
   #Y[Y==0] = NA
   #xs <- as(Y, "Incomplete")
   
   lam0 <- ifelse(is.na(lambda.init), lambda0.cov(Y, X) * lambda.factor, lambda.init) 
   #lam0 <- 40 
   lamseq <- seq(from=lam0, to=0, length=n.lambda)
   
   fits <- as.list(lamseq)
   ranks <- as.integer(lamseq)
   
   
   rank.max <- rank.init
   #warm <- NULL
   best_estimates = NA
   best_fit <- list(error=Inf, rank_A=NA, rank_B=NA, lambda=NA, rank.max=NA)
   counter <- 1
   X.X = t(X) %*% X
   if(n1n2 == 2){
      n1n2 = svd(X)$d[1]
   }else if(n1n2 == 3){
      n1n2 = nrow(Y) * ncol(Y)
   }
   if(lambda1 !=0) print("Setting lambda 1 to 0 to preseve the idempotent property of the Hat matrix")
   #H = solve(X.X +  diag(n1n2*lambda1, ncol(X))) %*% t(X)
   Q <- qr.Q(Matrix::qr(X)) #[,1:p]
   H <- Q %*% t(Q)
   
   
   for(i in seq(along=lamseq)) {
      fiti <- simpute.als.fit_splr(y = Y, yvalid = Y.valid, X = X, H = H,
                                   trace=F, J=rank.max, thresh=thresh, lambda=lamseq[i],
                                   warm.start = warm, patience=patience, maxit=100)
                                   
      
      # get test estimates and test error
      v=as.matrix(fiti$v)
      vd=v*outer(rep(1,nrow(v)),fiti$d)
      soft_estim = fiti$u %*% t(vd)  + X %*% fiti$beta.estim
      err = test_error(soft_estim[W==0], Y.valid)
      rank <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
      #----------------------------
      warm <- fiti # warm start for next 
      if(trace==TRUE)
         print(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
                       i, lamseq[i], rank.max, rank, err))
      #-------------------------
      # register best fir
      if(err < best_fit$error){
         best_fit$error = err
         best_fit$rank_A = qr(soft_estim)$rank
         best_fit$rank_B = rank
         best_fit$lambda = lamseq[i]
         best_fit$rank.max = rank.max
         best_estimates = soft_estim
         best_beta = fiti$beta.estim
         best_B = fiti$u %*% t(vd)
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
   if(print.best==TRUE) print(best_fit)
   best_fit$A_hat = best_estimates
   best_fit$B_hat = best_B
   best_fit$beta_hat = best_beta
   best_fit$last.fit = fiti
   return(best_fit)
}
#---