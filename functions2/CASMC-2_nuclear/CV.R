CASMC2_cv <-
   CASMC_Nuclear_cv <-
   function(y_train,
            # y_train is expected to be Incomplete
            X,
            y_valid,
            # y_valid is a vector
            W_valid,
            y = NULL,
            use_warmstart = TRUE,
            M_cv_param = list(
               rank.init = 2,
               rank.limit = 30,
               rank.step = 2,
               pct = 0.98,
               lambda.factor = 1 / 4,
               lambda.init = NULL,
               n.lambda = 20,
               early.stopping = 1
            ),
            beta_cv_param = list(
               rank.init = 2,
               rank.limit = qr(X)$rank,
               rank.step = 2,
               pct = 0.98,
               lambda.factor = 20,
               lambda.init = NULL,
               n.lambda = 20,
               early.stopping = 1
            ),
            error_function = utils$error_metric$rmse,
            
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
            print.best = TRUE,
            quiet = FALSE,
            track = FALSE,
            # initial values.
            warm = NULL,
            step3 = TRUE,
            seed = NULL) {
      if (!is.null(seed))
         set.seed(seed)
      
      #-----------------------------------------------------------------------
      # step 1: use highest rank of beta and no regularization to find optimal
      # parameters for M (CASMC2_cv_M)
      # step 2: Find parameters of beta without refitting M
      # step 3 [optional]: Try again with M to optimize.
      #------------------------------------------------------------------------
      
      # step 1:
      fiti = CASMC2_cv_M(
         y_train = y_train,
         X = X,
         y_valid = y_valid,
         W_valid = W_valid,
         y = NULL,
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
      )$fit
      M_param <- list(lambda = fiti$lambda.M, J = fiti$J)
      if(use_warmstart) warm = fiti
      #---------------------------------------------------------
      # step 2:
      fiti = CASMC2_cv_beta(
         y_train = y_train,
         X = X,
         y_valid = y_valid,
         W_valid = W_valid,
         y = NULL,
         J = M_param$J,
         lambda.M = M_param$lambda,
         warm = warm,
         thresh = thresh,
         maxit = maxit,
         # cv params
         lambda.factor = beta_cv_param$lambda.factor,
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
      beta_param <-
         list(lambda = fiti$lambda.beta, r = fiti$rank.max)
      if(use_warmstart) warm = fiti$fit
      #---------------------------------------------------------
      # step 3
      if (step3) {
         best_fit = CASMC2_cv_M(
            y_train = y_train,
            X = X,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            r = beta_param$r,
            lambda.beta = beta_param$lambda,
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
         )$fit
         M_param <- list(lambda = fiti$lambda.M, J = fiti$J)
         
      } else{
         best_fit <-
            CASMC2_fit(
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
               warm.start = warm,
               trace.it = trace,
               thresh = thresh,
               maxit = maxit
            )
         best_fit$beta = list(u = best_fit$ub,
                              v = best_fit$vb,
                              d = best_fit$db ^ 2)
      }
      #--------------------------------------------------------
      list(
         fit = best_fit,
         hparams = data.frame(
            lambda.beta = beta_param$lambda,
            rank.beta = beta_param$r,
            lambda.M = M_param$lambda,
            rank.M = M_param$J
         )
      )
   }
