CASMC_cv_holdout_v2 <-
   function(Y,
            X_r,
            Y_valid,
            W_valid,
            y = NULL,
            lambda.factor = 1 / 4,
            lambda.init = NA,
            n.lambda = 20,
            trace = FALSE,
            print.best = TRUE,
            tol = 1,
            thresh = 1e-6,
            maxit = 100,
            rank.init = 2,
            rank.limit = 30,
            rank.step = 2,
            warm = NULL,
            quiet = FALSE) {
      #' Input: Y training partially observed with missing set to 0
      #' X_r list
      #' Y_valid  vector of validation values
      #' W_valid matrix where W_valid = 0 for validation and 1 for train/test
      #' y (optional) whole training set including the validation for one last fit at the end.
      #' --------------------------------------------------------------------
      lam0 <-
         ifelse(is.na(lambda.init),
                lambda0.cov_splr(Y, X_r$svdH) * lambda.factor,
                lambda.init)
      lamseq <- seq(from = lam0,
                    to = 0,
                    length = n.lambda)
      #----------------------------------------------------
      stopifnot(inherits(Y, "dgCMatrix"))
      #Y[Y == 0] = NA
      #m = dim(Y)[2]
      #Y <- as(Y, "Incomplete")
      #valid_ind = W_valid == 0
      W_valid[W_valid == 1] = NA
      W_valid[W_valid == 0] =  1
      W_valid <- as(W_valid, "Incomplete")
      virow = W_valid@i
      vpcol = W_valid@p
      W_valid = NULL
      #-------------
      Xterms = GetXterms(X_r$X)
      #-----------------------------------------------------------------------
      rank.max <- rank.init
      best_fit <-
         list(
            error = Inf,
            rank_A = NA,
            rank_B = NA,
            lambda = NA,
            rank.max = NA
         )
      counter <- 0
      #---------------------------------------------------------------------
      for (i in seq(along = lamseq)) {
         fiti <-
            CASMC_fit_v3(
               y = Y,
               X = X_r$X,
               svdH = X_r$svdH,
               trace = F,
               J = rank.max,
               thresh = thresh,
               lambda = lamseq[i],
               init = "naive",
               final.svd = T,
               maxit = maxit,
               warm.start = warm,
               Xterms = Xterms
            )
         
         #--------------------------------------------------------------
         # predicting validation set and xbetas for next fit:
         
         Beta = fast.svd(as.matrix(fiti$beta))
         XbetaValid = suvC(X_r$X %*% Beta$v, t(Beta$d * t(Beta$u)), virow, vpcol)
         MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), virow, vpcol)
         #--------------------------------------------
         err = test_error(MValid + XbetaValid, Y_valid)
         rank <-
            sum(round(fiti$d, 4) > 0) # number of positive sing.values
         
         
         #---------------------------------------------------------------------
         #----------------------------
         warm <- fiti # warm start for next
         if (trace == TRUE)
            print(sprintf(
               paste0(
                  "%2d lambda=%9.5g, rank.max = %d  ==>",
                  " rank = %d, error = %.5f, niter/fit = %d"
               ),
               i,
               lamseq[i],
               rank.max,
               rank,
               err,
               fiti$n_iter
            ))
         #-------------------------
         # register best fir
         if (err < best_fit$error) {
            best_fit$error = err
            best_fit$rank_B = rank
            best_fit$lambda = lamseq[i]
            best_fit$rank.max = rank.max
            best_fit$fit = fiti
            best_fit$iter = i
            counter = 0
         } else
            counter = counter + 1
         if (counter >= tol) {
            if (trace)
               print(sprintf(
                  "Performance didn't improve for the last %d iterations.",
                  counter
               ))
            break
         }
         # compute rank.max for next iteration
         rank.max <- min(rank + rank.step, rank.limit)
      }
      # fit one last time full model, if the train/valid is provided
      if (!is.null(y)) {
         stopifnot(inherits(y, "dgCMatrix"))
         #y[y == 0] = NA
         #y <- as(y, "Incomplete")
         best_fit$fit <-
            CASMC_fit_v3(
               y = y,
               X = X_r$X,
               svdH = X_r$svdH,
               trace = F,
               J = best_fit$rank.max,
               thresh = thresh,
               lambda = best_fit$lambda,
               init = "naive",
               final.svd = T,
               maxit = maxit,
               warm.start = best_fit$fit1,
               Xterms = Xterms
            )
         
      }
      return(best_fit)
   }
#---
