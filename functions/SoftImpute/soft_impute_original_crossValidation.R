simpute.cv <- function(Y,
                         W,
                         A,
                         n.lambda = 20,
                         trace = FALSE,
                         print.best = TRUE,
                         tol = 5,
                         thresh = 1e-6,
                         rank.init = 10,
                         rank.limit = 50,
                         rank.step = 2,
                         maxit = 300,
                         biscale = FALSE) {
 # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
 Y[Y == 0] = NA
 #xs <- as(Y, "Incomplete")
 lam0 <- lambda0(Y)
 #lam0 <- 40
 lamseq <- seq(from = lam0,
               to = 0,
               length = n.lambda)
 
 fits <- as.list(lamseq)
 ranks <- as.integer(lamseq)
 
 
 rank.max <- rank.init
 warm <- NULL
 best_estimates = NA
 best_fit <-
  list(
   error = Inf,
   rank_A = NA,
   rank_B = NA,
   lambda = NA,
   rank.max = NA
  )
 counter <- 1
 if (biscale) {
  Y = biScale(Y, col.scale = F, row.scale = F)
 }
 for (i in seq(along = lamseq)) {
  fiti <-
   softImpute(
    Y,
    lambda = lamseq[i],
    rank.max = rank.max,
    warm = warm,
    thresh = thresh,
    maxit = maxit
   )
  
  # compute rank.max for next iteration
  rank <-
   sum(round(fiti$d, 4) > 0) # number of positive sing.values
  rank.max <- min(rank + rank.step, rank.limit)
  
  # get test estimates and test error
  #soft_estim = complete(Y, fiti)
  if (biscale) {
   soft_estim = complete(Y, fiti, unscale = TRUE)
  } else
   soft_estim = fiti$u %*% (fiti$d * t(fiti$v))
  err = test_error(soft_estim[W == 0], A[W == 0])
  #----------------------------
  warm <- fiti # warm start for next
  if (trace == TRUE)
   cat(
    sprintf(
     "%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
     i,
     lamseq[i],
     rank.max,
     rank,
     err
    )
   )
  #-------------------------
  # register best fir
  if (err < best_fit$error) {
   best_fit$error = err
   #best_fit$rank_A = qr(soft_estim)$rank
   best_fit$rank_B = rank
   best_fit$lambda = lamseq[i]
   best_fit$rank.max = rank.max
   #best_estimates = soft_estim
   counter = 1
  } else
   counter = counter + 1
  if (counter >= tol) {
   if (trace | print.best)
    cat(sprintf(
     "Performance didn't improve for the last %d iterations.",
     counter
    ))
   break
  }
 }
 #----------------------------------
 # one final fit on the whole data:
 if (biscale)
  A = biScale(A, col.scale = TRUE, row.scale = TRUE)
 fiti <-
  softImpute(
   A,
   lambda = best_fit$lambda,
   rank.max = best_fit$rank.max,
   warm = warm,
   thresh = thresh,
   maxit = maxit
  )
 if (biscale) {
  best_fit$A_hat = complete(A, fiti, unscale = TRUE)
 } else
  best_fit$A_hat =  fiti$u %*% (fiti$d * t(fiti$v))
 best_fit$rank_A = qr(best_fit$A_hat)$rank
 #----------------------------------------
 if (print.best == TRUE)
  print(best_fit)
 best_fit$B_hat = NA
 best_fit$beta_hat = NA
 return(best_fit)
}
