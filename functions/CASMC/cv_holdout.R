# CASMC_cv_holdout_with_r <-
#    function(y_train,
#             X_r,
#             y_valid,
#             W_valid,
#             r_min = 0,
#             y = NULL,
#             error_function = RMSE_error,
#             lambda.factor = 1 / 4,
#             lambda.init = NULL,
#             n.lambda = 20,
#             trace = FALSE,
#             print.best = TRUE,
#             early.stopping = 1,
#             thresh = 1e-6,
#             maxit = 100,
#             rank.init = 2,
#             rank.limit = 30,
#             rank.step = 2,
#             warm = NULL,
#             track_r = FALSE,
#             quiet = FALSE) {
#       r_seq <- (X_r$rank):(r_min)#(X_r$rank):(r_min)
#       Xterms = GetXterms(X_r$X)
#       best_score = Inf
#       best_fit = NULL
#       warm = NULL
#
#       for (r in r_seq) {
#          fiti <- CASMC_cv_holdout(
#             y_train = y_train,
#             X_r = X_r,
#             y_valid = y_valid,
#             W_valid = W_valid,
#             y = y,
#             Xterms = Xterms,
#             r = r,
#             error_function = error_function,
#             lambda.factor = lambda.factor,
#             lambda.init = lambda.init,
#             n.lambda = n.lambda,
#             trace = trace,
#             print.best = print.best,
#             early.stopping = early.stopping,
#             thresh = thresh,
#             maxit = maxit,
#             rank.init = rank.init,
#             rank.limit = rank.limit,
#             rank.step = rank.step,
#             warm = warm,
#             quiet = quiet
#          )
#          #warm = fiti$fit
#          if (fiti$error < best_score) {
#             best_score = fiti$error
#             best_fit = fiti
#          }
#          if (track_r)
#             print(paste(r, "-", fiti$error))
#       }
#       return(best_fit)
#
#    }

#------------------------------------------------------------------------------------------
CASMC_cv_holdout_with_r <-
   function(y_train,
            X_r,
            y_valid,
            W_valid,
            y = NULL,
            r_min = 0,
            r_max = X_r$rank,
            error_function = error_metric$rmse,
            lambda.factor = 1 / 4,
            lambda.init = NULL,
            n.lambda = 20,
            trace = FALSE,
            print.best = TRUE,
            early.stopping = 1,
            thresh = 1e-6,
            maxit = 100,
            rank.init = 2,
            rank.limit = 30,
            rank.step = 2,
            lambda.a = 0, 
            S.a = NULL,
            lambda.b = 0, 
            S.b=NULL,
            warm = NULL,
            track_r = FALSE,
            max_cores = 12,
            pct = 0.98,
            quiet = FALSE,
            seed = NULL) {
      r_seq <- (max(r_min, 0)):(min(X_r$rank, r_max))
      Xterms = GetXterms(X_r$X)
      best_score = Inf
      best_fit = NULL
      num_cores = min(max_cores, length(r_seq))
      print(paste("Running on", num_cores, "cores."))
      results <- mclapply(r_seq, function(r) {
         CASMC_cv_holdout(
            y_train = y_train,
            X_r = X_r,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            Xterms = Xterms,
            r = r,
            error_function = error_function,
            lambda.factor = lambda.factor,
            lambda.init = lambda.init,
            n.lambda = n.lambda,
            trace = trace,
            print.best = print.best,
            early.stopping = early.stopping,
            thresh = thresh,
            maxit = maxit,
            rank.init = rank.init,
            rank.limit = rank.limit,
            rank.step = rank.step,
            pct = pct,
            warm = NULL,
            quiet = quiet,
            seed = seed,
            lambda.a = lambda.a, S.a=S.a,lambda.b = lambda.b, S.b=S.b
         )
      }, mc.cores = num_cores)
      
      
      best_fit <-
         results[[which.min(sapply(results, function(x)
            x$error))]]
      
      if (track_r) {
         sapply(results, function(x)
            print(paste(x$r, "-", x$error)))
      }
      
      return(best_fit)
      
      
   }





#---------------------------------------------------------------------------------------------
#' Covariate-Adjusted-Sparse-Matrix-completion
#' Cross-Validation function with holdout method (training / validation)
#'
#' @param y_train A sparse matrix of class incomplete (training)
#' @param X_r A list containing SvdH and X (see fit function for more details)
#' @param y_valid A vector of the validation set.
#' @param W_valid An indicator matrix where 0 = validation and 1 = train/test
#' @param y (optional) A sparse matrix of class Incomplete of the train and validation combined. A final fit step will be
#'                            applied if y is provided.
#' @return A list of u,d,v of M, Beta, and a vector of the observed Xbeta
#' @examples
#'  CASMC_fit(y,X,J=5)
#' @export
#'
CASMC_cv_holdout <-
   function(y_train,
            X_r,
            y_valid,
            W_valid,
            r = NULL,
            y = NULL,
            Xterms = NULL,
            error_function = error_metric$rmse,
            lambda.factor = 1 / 4,
            lambda.init = NULL,
            n.lambda = 20,
            trace = FALSE,
            print.best = TRUE,
            early.stopping = 1,
            thresh = 1e-6,
            maxit = 100,
            rank.init = 2,
            rank.limit = 30,
            rank.step = 2,
            lambda.a = 0, 
            S.a = NULL,
            lambda.b = 0, 
            S.b=NULL,
            warm = NULL,
            pct = 0.98,
            quiet = FALSE,
            seed = NULL) {
      if (!is.null(seed))
         set.seed(seed)
      # prepare the sequence of lambda (nuclear regularization hyperparameter)
      if (is.null(lambda.init))
         lambda.init <-
            lambda0.cov_splr(y_train, X_r$svdH) * lambda.factor
      lamseq <- seq(from = lambda.init,
                    to = 0,
                    length = n.lambda)
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
      if (is.null(Xterms))
         Xterms = GetXterms(X_r$X)
      #-----------------------------------------------------------------------
      rank.max <- rank.init
      best_fit <- list(error = Inf, r = r)
      counter <- 0
      #---------------------------------------------------------------------
      for (i in seq(along = lamseq)) {
         fiti <-
            CASMC_fit(
               y = y_train,
               X = X_r$X,
               #svdH = X_r$svdH,
               Xterms = Xterms,
               r = r,
               J = rank.max,
               lambda = lamseq[i],
               warm.start = warm,
               trace.it = F,
               thresh = thresh,
               lambda.a = lambda.a, 
               S.a = S.a,
               lambda.b = lambda.b, 
               S.b=S.b,
               init = "naive",
               final.svd = T,
               maxit = maxit
            )
         
         #--------------------------------------------------------------
         # predicting validation set and xbetas for next fit:
         Beta = fiti$Beta
         XbetaValid = suvC(X_r$X %*% Beta$v, t(Beta$d * t(Beta$u)), virow, vpcol)
         MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), virow, vpcol)
         #--------------------------------------------
         err = error_function(MValid + XbetaValid, y_valid)
         rank <- sum(round(fiti$d, 4) > 0)
         # newly added, to be removed later
         var_explained = fiti$d ^ 2 / sum(fiti$d ^ 2)
         cum_var = cumsum(var_explained)
         rank <- rank2 <- which(cum_var >= pct)[1]
         #print( fiti$d)
         #rank <- rank2 <- min(rank, 2, na.rm=TRUE)
         warm <- fiti # warm start for next
         #print(paste(rank,"-",rank2))
         #---------------------------------------------------------------------
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
            best_fit$rank_M = rank
            best_fit$lambda = lamseq[i]
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
         rank.max <- min(rank2 + rank.step, rank.limit)
      }
      # fit one last time full model, if the train/valid is provided
      if (!is.null(y)) {
         stopifnot(inherits(y, "dgCMatrix"))
         best_fit$fit <-
            CASMC_fit(
               y = y,
               X = X_r$X,
               #svdH = X_r$svdH,
               Xterms = Xterms,
               r = r,
               J = best_fit$rank.max,
               lambda = best_fit$lambda,
               warm.start = best_fit$fit1,
               thresh = thresh,
               maxit = maxit,
               trace.it = F,
               lambda.a = lambda.a, 
               S.a = S.a,
               lambda.b = lambda.b, 
               S.b=S.b,
               final.svd = T
            )
         
      }
      return(best_fit)
   }
#---
