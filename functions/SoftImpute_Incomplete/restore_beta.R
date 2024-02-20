require(svd)
restore.beta <- function(Y, fits, thresh=1e-6, maxit=300, trace.it=T){
   fits$M = fits$u %*% (fits$d * t(fits$v))
   Y = Y_train
   ynas = Y == 0
   Mmis = fits$M[ynas]
   Y[ynas] = Mmis
   Y[!ynas] = Y[!ynas] - fits$xbeta.obs
   XY = MASS::ginv(X_r$X) %*% Y
   rank = fits$rank
   
   old_beta = NULL
   ratio = Inf
   iter = 1
   obj = NA
   
   while((ratio > thresh)&(iter < maxit)){
      
      svd_beta = fast.svd(XY)     #propack.svd(Y, neig=rank)
      beta = svd_beta$u %*% (svd_beta$d * t(svd_beta$v))
      
      Y[ynas] = Mmis - (X_r$X %*% beta)[ynas]
      XY = MASS::ginv(X_r$X) %*% Y
      
      
      if(iter > 1) ratio = sum( (beta - old_beta)^2 ) / sum((old_beta)^2)
      if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio,"\n")
      iter = iter + 1
      old_beta = beta
   }
   
   
   fits$A= Y
}