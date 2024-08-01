simpute.cv <- function(Y_train,
                       y_valid,
                       W_valid,
                       y,
                       n.lambda = 20,
                       trace = FALSE,
                       print.best = TRUE,
                       tol = 5,
                       thresh = 1e-6,
                       rank.init = 10,
                       rank.limit = 50,
                       rank.step = 2,
                       maxit = 300,
                       test_error = utils$error_metric$rmse,
                       seed = NULL) {
  # W: validation only wij=0. For train and test make wij=1. make Yij=0 for validation and test. Aij=0 for test only.
  if(!is.null(seed))
    set.seed(seed)
  valid_ind <- W_valid == 0
  y[y == 0] = NA
  Y_train[Y_train == 0] = NA
  lam0 <- softImpute::lambda0(y)
  lamseq <- seq(from = lam0,
                to = 0,
                length = n.lambda)
  
  rank.max <- rank.init
  warm <- NULL
  best_fit <-
    list(
      error = Inf,
      rank_M = NA,
      lambda = NA,
      rank.max = NA
    )
  counter <- 1
  
  for (i in seq(along = lamseq)) {
    fiti <-
      softImpute::softImpute(
        Y_train,
        lambda = lamseq[i],
        rank.max = rank.max,
        warm.start = warm,
        thresh = thresh,
        maxit = maxit
      )
    
    # compute rank.max for next iteration
    rank <-
      sum(round(fiti$d, 4) > 0) # number of positive sing.values
    rank.max <- min(rank + rank.step, rank.limit)
    
    
    soft_estim = fiti$u %*% (fiti$d * t(fiti$v))
    err = test_error(soft_estim[valid_ind], y_valid)
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
      best_fit$rank_M = rank
      best_fit$lambda = lamseq[i]
      best_fit$rank.max = rank.max
      counter = 1
    } else
      counter = counter + 1
    if (counter >= tol) {
      if (trace || print.best)
        cat(sprintf(
          "Performance didn't improve for the last %d iterations.",
          counter
        ))
      break
    }
  }
  #----------------------------------------
  if (print.best == TRUE)
    print(best_fit)
  #----------------------------------
  # one final fit on the whole data:
  

    fiti <-
    softImpute::softImpute(
      y,
      lambda = best_fit$lambda,
      rank.max = best_fit$rank.max,
      warm.start = warm,
      thresh = thresh,
      maxit = maxit
    )
  
  best_fit$estimates =  fiti$u %*% (fiti$d * t(fiti$v))
  best_fit$rank_M = qr(best_fit$estimates)$rank
  #---------------------------------------------------------------
  return(best_fit)
}
