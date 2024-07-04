CASMC2_cv_kf <-
   function(Y,
            # y_train is expected to be Incomplete
            X,
            y,
            obs_mask,
            n_folds = 5,
            M_cv_param = list(
               rank.init = 2,
               rank.limit = 30,
               rank.step = 2,
               pct = 0.98,
               lambda.factor = 1/4,
               lambda.init = NULL,
               n.lambda = 20, 
               early.stopping = 1
            ),
            beta_cv_param = list(
               rank.init = 2,
               rank.limit = qr(X)$rank,
               rank.step = 2,
               pct = 0.98,
               lambda.multi.factor = 20,
               lambda.init = NULL,
               n.lambda = 20, 
               early.stopping = 1
            ),
            
            # y: a final full-fit if provided. Expected to be Incomplete
            error_function = error_metric$rmse,
            # laplacian parameters
            lambda.a = 0,
            S.a = NULL,
            lambda.b = 0,
            S.b = NULL,
            # stopping criteria
            thresh = 1e-6,
            maxit = 100,
            # trace parameters
            trace = FALSE,
            track = FALSE,
            print.best = TRUE,
            quiet = FALSE,
            # initial values.
            warm = NULL,
            max_cores = 12,
            # seed
            seed = NULL)
   {
      if (!is.null(seed))
         set.seed(seed)
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
         mask_tmp[valid_ind] = 0
         mask_tmp[!valid_ind] = 1
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
      num_cores = n_folds
      if (num_cores > max_cores)
         num_cores <-
         min(max_cores, ceiling(n_folds / 2))
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
      print("starting 1")
      fits <-
         mclapply(1:n_folds, function(fold) {
            fdat = fold_data[[fold]]
            
            
            fiti = CASMC2_cv_M(
               y_train = fdat$Y_train,
               X = X,
               y_valid = fdat$Y_valid,
               W_valid = fdat$valid_mask,
               r = beta_cv_param$rank.limit,
               lambda.beta = 0,
               warm = warm,
               thresh = thresh,
               maxit = maxit,
               # cv params
               lambda.factor = M_cv_param$lambda.factor,
               lambda.init = M_cv_param$lambda.init,
               n.lambda = M_cv_param$n.lambda,
               rank.init = M_cv_param$rank.init,
               rank.limit = M_cv_param$rank.limit,
               rank.step = M_cv_param$rank.step,
               pct = M_cv_param$pct,
               early.stopping = M_cv_param$early.stopping,
               # error and others
               error_function = error_function,
               trace = trace,
               print.best = print.best,
               quiet = quiet,
               lambda.a = lambda.a,
               S.a = S.a,
               lambda.b = lambda.b,
               S.b = S.b,
               seed = NULL
            )
            
            list(
               fit = fiti$fit,
               error = fiti$error,
               lambda.M = fiti$fit$lambda.M,
               J = fiti$fit$J
            )
            
         }, mc.cores =  num_cores, mc.preschedule = TRUE)
      #------------------------------------
      print("end of 1")
      # Post-step 1: Average out the M parameters
      lambda.M_values <- sapply(fits, function(s) s$lambda.M)
      J_values <- sapply(fits, function(s) s$J)
      M_param <- list(lambda = mean(lambda.M_values), J = round(max(J_values)))
      #---------------------------------------------------------
      # step 2:
      print("startin 2")
      fits <- mclapply(1:n_folds, function(fold) {
         fdat = fold_data[[fold]]
         warm = fits[[fold]]$fit
         
         fiti = CASMC2_cv_beta(
            y_train = fdat$Y_train,
            X = X,
            y_valid = fdat$Y_valid,
            W_valid = fdat$valid_mask,
            J = M_param$J,
            lambda.M = M_param$lambda,
            warm = warm,
            thresh = thresh,
            maxit = maxit,
            # cv params
            lambda.multi.factor = beta_cv_param$lambda.multi.factor,
            lambda.init = beta_cv_param$lambda.init,
            n.lambda = beta_cv_param$n.lambda,
            rank.init = beta_cv_param$rank.init,
            rank.limit = beta_cv_param$rank.limit,
            rank.step = beta_cv_param$rank.step,
            pct = beta_cv_param$pct,
            early.stopping = beta_cv_param$early.stopping,
            # error and others
            error_function = error_function,
            trace = trace,
            print.best = print.best,
            quiet = quiet,
            lambda.a = lambda.a,
            S.a = S.a,
            lambda.b = lambda.b,
            S.b = S.b,
            seed = NULL
         )
         
         list(
            fit = fiti$fit,
            error = fiti$error,
            lambda.beta = fiti$lambda.beta,
            r = fiti$rank.max
         )
         
      },  mc.cores =  num_cores, mc.preschedule = TRUE)
      #------------------------------------
      print("end of 2")
      # Post-step 2: Average out the beta parameters
      lambda.beta_values <- sapply(fits, function(s) s$lambda.beta)
      r_values <- sapply(fits, function(s) s$r)
      beta_param <- list(lambda = mean(lambda.beta_values), r = round(max(r_values)))
      #---------------------------------------------------------
      # step 3:
      print("3 starting")
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
               lambda.beta = beta_param$r,
               warm = warm,
               thresh = thresh,
               maxit = maxit,
               # cv params
               lambda.factor = M_cv_param$lambda.factor,
               lambda.init = M_cv_param$lambda.init,
               n.lambda = M_cv_param$n.lambda,
               rank.init = M_cv_param$rank.init,
               rank.limit = M_cv_param$rank.limit,
               rank.step = M_cv_param$rank.step,
               pct = M_cv_param$pct,
               early.stopping = M_cv_param$early.stopping,
               # error and others
               error_function = error_function,
               trace = trace,
               print.best = print.best,
               quiet = quiet,
               lambda.a = lambda.a,
               S.a = S.a,
               lambda.b = lambda.b,
               S.b = S.b,
               seed = NULL
            )
            
            list(
               fit = fiti$fit,
               error = fiti$error,
               lambda.M = fiti$fit$lambda.M,
               J = fiti$fit$J
            )
            
         }, mc.cores =  num_cores, mc.preschedule = TRUE)
      #------------------------------------
      print("end of 3")
      # post-step 3: average out the M parameters
      lambda.M_values <- sapply(fits, function(s) s$lambda.M)
      J_values <- sapply(fits, function(s) s$J)
      M_param <- list(lambda = mean(lambda.M_values), J = round(max(J_values)))
      #-----------------------------------------------------
      #output
      
      
      M <- matrix(0, nrow(Y), ncol(Y))
      beta <- matrix(0, ncol(X), ncol(Y))
      for (fold in 1:n_folds) {
         M <- M + unsvd(fits[[fold]]$fit)
         beta <- beta + unsvd(fits[[fold]]$fit$beta)
      }
      M <- svd_trunc_simple(M / n_folds, M_param$J)
      beta <- svd_trunc_simple(beta / n_folds, beta_param$r)
      warm <- list(
         u = M$u,
         d = M$d,
         v = M$v,
         
         ub = beta$u,
         db = sqrt(beta$d),
         vb = beta$v
      )
      
      fit <- CASMC2_fit(
         y = y,
         X = X,
         J = M_param$J,
         lambda.M = M_param$lambda,
         r = beta_param$r,
         lambda.beta = beta_param$lambda,
         lambda.a = lambda.a,
         S.a = S.a,
         lambda.b = lambda.b,
         S.b = S.b,
         warm.start = warm,#fits[[1]]$fit,
         trace.it = T,
         thresh = thresh,
         maxit = maxit
      )
      fit$beta <- list(
         u = fit$ub,
         d = fit$db^2,
         v = fit$vb
      )
      
      
      #--------------------------------------------------------
      list(
         fits = fits,
         M = unsvd(fit),
         beta = unsvd(fit$beta),
         M_param = M_param,
         beta_param = beta_param
      )
      
   }
