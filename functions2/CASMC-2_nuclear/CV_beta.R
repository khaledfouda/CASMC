CASMC2_cv_beta <-
 function(y_train,
          # y_train is expected to be Incomplete
          X,
          y_valid,
          # y_valid is a vector
          W_valid,
          J,
          lambda.M,
          y = NULL,
          # y: a final full-fit if provided. Expected to be Incomplete
          error_function = utils$error_metric$rmse,
          # tuning parameters for lambda
          lambda.factor = 20,
          lambda.init = NULL,
          n.lambda = 20,
          # tuning parameters for J
          rank.init = 1,
          rank.limit = qr(X)$rank,
          rank.step = 1,
          pct = 0.98,
          # laplacian parameters
          lambda.a = 0,
          S.a = NULL,
          lambda.b = 0,
          S.b = NULL,
          # stopping criteria
          early.stopping = 1,
          thresh = 1e-6,
          maxit = 100,
          # trace parameters
          trace = FALSE,
          print.best = TRUE,
          quiet = TRUE,
          # initial values.
          warm = NULL,
          # seed
          seed = NULL) {
  if (!is.null(seed))
   set.seed(seed)
  
  # prepare the sequence of lambda (nuclear regularization hyperparameter)
  if (is.null(lambda.init)) {
   # lamseq = sqrt((ncol(y_train) * ncol(X)) / (nrow(y_train))) *
   #  seq(lambda.factor, .Machine$double.eps, length.out = n.lambda)
   
   lambda.init <-
     utils$lambda.beta.max(y_train, utils$reduced_hat_decomp.H(X)) * lambda.factor
   lamseq <- seq(from = lambda.init,
                 to = .Machine$double.eps,
                 length = n.lambda)
   
   
  } else{
   lamseq = seq(lambda.init, .Machine$double.eps, length.out = n.lambda)
  }
  #----------------------------------------------------
  stopifnot(inherits(y_train, "dgCMatrix"))
  # we only need the indices for validation from W_valid
  W_valid[W_valid == 1] = NA
  W_valid[W_valid == 0] =  1
  W_valid <- as(W_valid, "Incomplete")
  virow = W_valid@i
  vpcol = W_valid@p
  W_valid = NULL
  #------------------------------------------------
  #-----------------------------------------------------------------------
  rank.max <- rank.init
  best_fit <- list(error = Inf)
  counter <- 0
  #---------------------------------------------------------------------
  for (i in seq(along = lamseq)) {
   fiti <-
    CASMC2_fit(
     y = y_train,
     X = X,
     J = J,
     r = rank.max,
     lambda.M = lambda.M,
     lambda.beta = lamseq[i],
     lambda.a = lambda.a,
     S.a = S.a,
     lambda.b = lambda.b,
     S.b = S.b,
     warm.start = warm,
     trace.it = F,
     thresh = thresh,
     maxit = maxit
    )
   
   #--------------------------------------------------------------
   # predicting validation set and xbetas for next fit:
   XbetaValid = suvC(as.matrix(X %*% fiti$ub),
                     as.matrix(utils$UD(fiti$vb, fiti$db ^ 2)),
                     virow,
                     vpcol)
   MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), virow, vpcol)
   #--------------------------------------------
   err = error_function(MValid + XbetaValid, y_valid)
   #rank <- sum(round(fiti$d, 4) > 0)
   # newly added, to be removed later
   var_explained = fiti$db ^ 2 / sum(fiti$db ^ 2)
   cum_var = cumsum(var_explained)
   rank  <- which(cum_var >= pct)[1]
   
   warm <- fiti # warm start for next
   #---------------------------------------------------------------------
   if (trace == TRUE)
    print(sprintf(
     paste0(
      "[Beta]: %2d lambda=%9.5g, rank.max = %d  ==>",
      " rank = %d, error = %.5f, niter/fit = %d [M(J=%d, lambda=%.3f)]"
     ),
     i,
     lamseq[i],
     rank.max,
     rank,
     err,
     fiti$n_iter,
     J,
     lambda.M
    ))
   #-------------------------
   # register best fir
   if (err < best_fit$error) {
    best_fit$error = err
    best_fit$rank.beta = rank
    best_fit$lambda.beta = lamseq[i]
    best_fit$rank.max = rank.max
    best_fit$fit = fiti
    best_fit$iter = i
    counter = 0
   } else
    counter = counter + 1
   if (counter >= early.stopping) {
    if (trace)
     print(
      sprintf(
       "Early stopping. Reached Peak point. Performance didn't improve for the last %d iterations.",
       counter
      )
     )
    break
   }
   # compute rank.max for next iteration
   rank.max <- min(rank + rank.step, rank.limit)
  }
  # fit one last time full model, if the train/valid is provided
  if (!is.null(y)) {
   stopifnot(inherits(y, "dgCMatrix"))
   best_fit$fit <-
    CASMC2_fit(
     y = y,
     X = X,
     r = best_fit$rank.max,
     lambda.beta = best_fit$lambda.beta,
     J = J,
     lambda.M = lambda.M,
     lambda.a = lambda.a,
     S.a = S.a,
     lambda.b = lambda.b,
     S.b = S.b,
     warm.start = best_fit$fit,
     trace.it = F,
     thresh = thresh,
     maxit = maxit
    )
   
  }
  best_fit$lambda.M = lambda.M
  best_fit$J = J
  best_fit$lambda.a = lambda.a
  best_fit$lambda.b = lambda.b
  best_fit$fit$beta <- list(u = best_fit$fit$ub,
                            d = best_fit$fit$db ^ 2,
                            v = best_fit$fit$vb)
  return(best_fit)
 }
#---