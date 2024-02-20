require(svd)
restore.beta <- function(Y, M,  xbeta.S, rank, thresh=1e-6, maxit=300, trace.it=T){
   
   #fits$M = fits$u %*% (fits$d * t(fits$v))
   #Y = Y_train
   #xbeta.S = fits$xbeta.sparse
   
   ynas = Y == 0
   Mmis = M[ynas]
   Y[ynas] = Mmis
   Y = Y - xbeta.S
   #XY = MASS::ginv(X_r$X) %*% Y

   old_beta = beta = NULL
   ratio = Inf
   iter = 0
   obj = NA
   
   time1 = time2 = time3 = 0
   
   while((ratio > thresh)&(iter < maxit)){
      iter = iter + 1
      old_beta = beta
      start_time <- Sys.time()
      svd_beta = propack.svd(xbeta.S, neig=rank)
      time1 = time1 + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
      
      start_time <- Sys.time()
      beta = svd_beta$u %*% (svd_beta$d * t(svd_beta$v))
      xbeta.S[ynas] = beta[ynas]
      time2 = time2 + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
      
      start_time <- Sys.time()
      if(iter > 1) ratio = sum( (beta - old_beta)^2 ) / sum((old_beta)^2)
      time3 = time3 + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
      if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio,"\n")
      
   }
   print(paste(time1, time2, time3))
   list(xbeta = xbeta.S, A = Y + xbeta.S)
}
