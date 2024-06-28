CASMC2_cv2 <-
   function(Y,
            # y_train is expected to be Incomplete
            X,
            obs_mask,
            n_folds = 5,
            rank.beta.init = 1,
            rank.beta.step = 1,
            rank.beta.limit = qr(X)$rank,
            lambda.beta.length = 20,
            lambda.beta.max = NULL,
            # y: a final full-fit if provided. Expected to be Incomplete
            error_function = error_metric$rmse,
            # tuning parameters for lambda
            lambda.factor = 1 / 4,
            lambda.init = NULL,
            n.lambda = 20,
            # tuning parameters for J
            rank.M.init = 2,
            rank.M.limit = 30,
            rank.step = 2,
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
            quiet = FALSE,
            # initial values.
            warm = NULL,
            step3 = TRUE,
            # L2 parameters
            lambda.beta.grid = "default1",
            track = FALSE,
            max_cores = 12,
            # seed
            seed = NULL) {
      if (!is.null(seed))
         set.seed(seed)
      
      if (is.null(lambda.beta.max)) {
         if (identical(lambda.beta.grid, "default1"))
            lambda.beta.grid = sqrt((ncol(Y) * ncol(X)) / (nrow(Y))) *
               seq(20, .Machine$double.eps, length.out = lambda.beta.length)
         
         
         if (identical(lambda.beta.grid, "default2"))
            lambda.beta.grid = seq(propack.svd(naive_fit(Y, X, TRUE), 1)$d / 4,
                                   .Machine$double.eps,
                                   length.out = lambda.beta.length)
      } else
         lambda.beta.grid = seq(lambda.beta.max, 0, length.out = lambda.beta.length)
      
      #-----------------------------------------------------------------------
      # prepare the folds
      folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, obs_mask)
      fold_data <- lapply(1:n_folds, function(i) {
         valid_mask = folds[[i]]
         valid_ind = valid_mask == 0 & obs_mask == 1
         
         Y_train = Y * valid_mask
         Y_train[valid_mask == 0] = NA
         Y_train = as(Y_train, "Incomplete")
         
         Y_valid = Y[valid_ind]
         
         mask_tmp <- valid_mask
         valid_mask[valid_ind] = 1
         valid_mask[!valid_ind] = NA
         valid_mask <- as(valid_mask, "Incomplete")
         virow = valid_mask@i
         vpcol = valid_mask@p
         valid_mask <- NULL
         
         list(
            Y_train = Y_train,
            Y_valid = Y_valid,
            virow = virow,
            vpcol = vpcol,
            valid_mask = mask_tmp
         )
      })
      #------------------------------------------------
      
      
      num_cores = length(n_folds)
      if (num_cores > max_cores)
         num_cores <-
         min(max_cores, ceiling(length(n_folds) / 2))
      print(paste("Running on", num_cores, "cores."))
      
      #-----------------------------------------------------------------------
      # step 1: use highest rank of beta and no regularization to find optimal
      # parameters for M (CASMC2_cv_M)
      # step 2: Find parameters of beta without refitting M
      # step 3 [optional]: Try again with M to optimize.
      #------------------------------------------------------------------------
      
      # step 1:
      warm = NULL
      fits <- vector("list", n_folds)
      
      fits <-
         mclapply(1:n_folds, function(fold) {
            fdat = fold_data[[fold]]
            
            
            fiti = CASMC2_cv_M(
               y_train = fdat$Y_train,
               X = X,
               y_valid = fdat$Y_valid,
               W_valid = fdat$valid_mask,
               r = rank.beta.limit,
               lambda.beta = 0,
               
               error_function = error_function,
               lambda.factor = lambda.factor,
               lambda.init = lambda.init,
               n.lambda = n.lambda,
               trace = trace,
               print.best = print.best,
               thresh = thresh,
               maxit = maxit,
               rank.init = rank.M.init,
               rank.limit = rank.M.limit,
               rank.step = rank.step,
               pct = pct,
               warm = warm,
               lambda.a = lambda.a,
               S.a = S.a,
               lambda.b = lambda.b,
               S.b = S.b,
               quiet = quiet,
               seed = seed
            )
            
            list(
               fit = fiti$fit,
               error = fiti$error,
               lambda.M = fiti$fit$lambda.M,
               J = fiti$fit$J
            )
            
         }, mc.cores =  num_cores)
      #------------------------------------
      # post-step 1: average out the M parameters
      M_param <- list(lambda = mean(sapply(fits, function(s)
         s$lambda.M)),
         J = mean(sapply(fits, function(s)
            s$J)))
      #---------------------------------------------------------
      # step 2:
      fits <- mclapply(1:n_folds, function(fold) {
         fdat = fold_data[[fold]]
         warm = fits[[fold]]$fit
         
         rank.max = rank.beta.init
         counter = 0
         best_fit <- list(error = Inf)
         for (i in seq(along = lambda.beta.grid)) {
            fiti <- CASMC2_fit(
               y = fdat$Y_train,
               X = X,
               J = M_param$J,
               lambda.M = M_param$lambda,
               r = rank.max,
               lambda.beta = lambda.beta.grid[i],
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
                              as.matrix(UD(fiti$vb, fiti$db ^ 2)),
                              fdat$virow,
                              fdat$vpcol)
            MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), fdat$virow, fdat$vpcol)
            #--------------------------------------------
            err = error_function(MValid + XbetaValid, fdat$y_valid)
            # rank <- sum(round(fiti$db, 6) > 0)
            var_explained = fiti$db ^ 2 / sum(fiti$db ^ 2)
            cum_var = cumsum(var_explained)
            rank  <- which(cum_var >= pct)[1]
            warm <- fiti
            #-----------------------------------------
            if (track == TRUE)
               print(
                  sprintf(
                     paste0(
                        "%2d lambda.beta=%9.5g, rank.max.beta = %d  ==>",
                        " rank.max = %d, error = %.5f, niter/fit = %d [M(J=%d, lambda=%.3f)]"
                     ),
                     i,
                     lambda.beta.grid[i],
                     rank.max,
                     rank,
                     err,
                     fiti$n_iter,
                     M_param$J,
                     M_param$lambda
                  )
               )
            #-------------------------
            # register best fit
            if (err < best_fit$error) {
               best_fit$error = err
               best_fit$rank_beta = rank
               best_fit$lambda.beta = lambda.beta.grid[i]
               best_fit$rank.beta.max = rank.max
               best_fit$fit = fiti
               best_fit$iter = i
               counter = 0
            } else
               counter = counter + 1
            if (counter >= early.stopping) {
               if (track)
                  print(
                     sprintf(
                        "Early stopping. Reached Peak point. Performance didn't improve for the last %d iterations.",
                        counter
                     )
                  )
               break
            }
            # compute rank.max for next iteration
            rank.max <- min(rank + rank.beta.step, rank.beta.limit)
         }
         list(fit = best_fit$fit,
              error = best_fit$error,
              lambda.beta = best_fit$lambda.beta,
              r = best_fit$rank.beta.max)
         
      })
      #------------------------------------
      # post-step 2: average out the beta parameters
      beta_param <- list(lambda = mean(sapply(fits, function(s)
         s$lambda.beta)),
         r = mean(sapply(fits, function(s)
            s$r)))
      #---------------------------------------------------------
      # step 3:
      fits <-
         mclapply(1:n_folds, function(fold) {
            fdat = fold_data[[fold]]
            warm = fits[[fold]]$fit   
            
            fiti = CASMC2_cv_M(
               y_train = fdat$Y_train,
               X = X,
               y_valid = fdat$Y_valid,
               W_valid = fdat$valid_mask,
               y = y, 
               r = beta_param$r,
               lambda.beta = beta_param$lambda,
               
               error_function = error_function,
               lambda.factor = lambda.factor,
               lambda.init = lambda.init,
               n.lambda = n.lambda,
               trace = trace,
               print.best = print.best,
               thresh = thresh,
               maxit = maxit,
               rank.init = rank.M.init,
               rank.limit = rank.M.limit,
               rank.step = rank.step,
               pct = pct,
               warm = warm,
               lambda.a = lambda.a,
               S.a = S.a,
               lambda.b = lambda.b,
               S.b = S.b,
               quiet = quiet,
               seed = seed
            )
            
            list(
               fit = fiti$fit,
               error = fiti$error,
               lambda.M = fiti$fit$lambda.M,
               J = fiti$fit$J
            )
            
         }, mc.cores =  num_cores)
      #------------------------------------
      # post-step 3: average out the M parameters
      M_param <- list(lambda = mean(sapply(fits, function(s)
         s$lambda.M)),
         J = mean(sapply(fits, function(s)
            s$J)))
      M <- matrix(0, nrow(Y), ncol(Y))
      beta <- matrix(0, ncol(X), ncol(Y))
      for(fold in 1:n_folds){
         M <- M + unsvd(fits[[fold]]$fit)
      }
      
      #---------------------------------------------------------
      # step 3
      if (step3) {
         best_fit$fit = CASMC2_cv_M(
            y_train = y_train,
            X = X,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            r = best_fit$rank_beta,
            lambda.beta = best_fit$lambda.beta,
            
            error_function = error_function,
            lambda.factor = lambda.factor,
            lambda.init = lambda.init,
            n.lambda = n.lambda,
            trace = trace,
            print.best = print.best,
            #early.stopping = early.stopping,
            thresh = thresh,
            maxit = maxit,
            rank.init = rank.M.init,
            rank.limit = rank.M.limit,
            rank.step = rank.step,
            pct = pct,
            warm = best_fit$fit,
            lambda.a = lambda.a,
            S.a = S.a,
            lambda.b = lambda.b,
            S.b = S.b,
            quiet = quiet,
            seed = seed
         )$fit
      } else{
         best_fit$fit <-
            CASMC2_fit(
               y = y,
               X = X,
               J = M_param$J,
               lambda.M = M_param$lambda,
               r = best_fit$rank_beta,
               lambda.beta = best_fit$lambda.beta,
               lambda.a = lambda.a,
               S.a = S.a,
               lambda.b = lambda.b,
               S.b = S.b,
               warm.start = best_fit$fit,
               trace.it = T,
               thresh = thresh,
               maxit = maxit
            )
      }
      #--------------------------------------------------------
      
      best_fit$hparams = data.frame(
         lambda.beta = best_fit$lambda.beta,
         lambda.M = best_fit$fit$lambda.M,
         rank.beta = best_fit$rank_beta,
         rank.M = best_fit$fit$J
      )
      best_fit$fit$beta = list(u = best_fit$fit$ub,
                               v = best_fit$fit$vb,
                               d = best_fit$fit$db ^ 2)
      
      return(best_fit)
      
      
   }



# rank.max = rank.beta.init
# counter = 0
# best_fit <- list(error = Inf)
#
# for (i in seq(along = lambda.beta.grid)) {
#    fiti = CASMC2_cv_M(
#       y_train = y_train,
#       X = X,
#       y_valid = y_valid,
#       W_valid = W_valid,
#       y = y,
#       r = rank.max,
#       lambda.beta = lambda.beta.grid[i],
#
#       error_function = error_function,
#       lambda.factor = lambda.factor,
#       lambda.init = lambda.init,
#       n.lambda = n.lambda,
#       trace = trace,
#       print.best = print.best,
#       #early.stopping = early.stopping,
#       thresh = thresh,
#       maxit = maxit,
#       rank.init = rank.M.init,
#       rank.limit = rank.M.limit,
#       rank.step = rank.step,
#       pct = pct,
#       warm = warm,
#       #warm,
#       lambda.a = lambda.a,
#       S.a = S.a,
#       lambda.b = lambda.b,
#       S.b = S.b,
#       quiet = quiet,
#       seed = seed
#    )
#    err = fiti$error
#    fiti = fiti$fit
#    rank <- sum(round(fiti$db, 6) > 0)
#    # var_explained = diag(fiti$db) ^ 2 / sum(diag(fiti$db) ^ 2)
#    # cum_var = cumsum(var_explained)
#    # rank  <- which(cum_var >= pct)[1]
#
#    warm <- fiti # warm start for next
#
#    if (track == TRUE)
#       print(sprintf(
#          paste0(
#             "%2d lambda.beta=%9.5g, rank.max.beta = %d  ==>",
#             " rank.max = %d, error = %.5f, niter/fit = %d"
#          ),
#          i,
#          lambda.beta.grid[i],
#          rank.max,
#          rank,
#          err,
#          fiti$n_iter
#       ))
#    #-------------------------
#    # register best fit
#    if (err < best_fit$error) {
#       best_fit$error = err
#       best_fit$rank_beta = rank
#       best_fit$lambda.beta = lambda.beta.grid[i]
#       best_fit$rank.beta.max = rank.max
#       best_fit$fit = fiti
#       best_fit$iter = i
#       counter = 0
#    } else
#       counter = counter + 1
#    if (counter >= early.stopping) {
#       if (track)
#          print(
#             sprintf(
#                "Early stopping. Reached Peak point. Performance didn't improve for the last %d iterations.",
#                counter
#             )
#          )
#       break
#    }
#    # compute rank.max for next iteration
#    rank.max <- min(rank + rank.beta.step, rank.beta.limit)
#
# }
