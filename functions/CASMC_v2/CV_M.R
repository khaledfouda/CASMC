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
CASMC2_cv_M <-
 function(y_train,
          # y_train is expected to be Incomplete
          X,
          y_valid,
          # y_valid is a vector
          W_valid,
          y = NULL,
          # y: a final full-fit if provided. Expected to be Incomplete
          r = NULL,
          lambda.beta = 0,
          # provide this if you need rank restriction. if not null, L2 reg will be ignored
          svdH = NULL,
          # provide this if you need L2 regularization.
          error_function = error_metric$rmse,
          # tuning parameters for lambda
          lambda.factor = 1 / 4,
          lambda.init = NULL,
          n.lambda = 20,
          # tuning parameters for J
          rank.init = 2,
          rank.limit = 30,
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
          # seed
          seed = NULL) {
  if (!is.null(seed))
   set.seed(seed)
  if(is.null(svdH)) 
   svdH <-  reduced_hat_decomp.H(X)
  # prepare the sequence of lambda (nuclear regularization hyperparameter)
  if (is.null(lambda.init))
   lambda.init <-
    lambda0.cov_splr(y_train, svdH) * lambda.factor
  lamseq <- seq(from = lambda.init,
                to = 0,
                length = n.lambda)
  #---------------------------------------------------
  # if (!is.null(r) || is.null(Xterms)) {
  #  fit_function = CASMC_fit_rank
  # } else
  #  fit_function = CASMC_fit_L2
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
  best_fit <- list(error = Inf, r = r)
  counter <- 0
  #---------------------------------------------------------------------
  for (i in seq(along = lamseq)) {
   fiti <-
    CASMC2_fit2(
     y = y_train,
     X = X,
     J = rank.max,
     lambda.M = lamseq[i],
     lambda.beta = lambda.beta,
     svdH = NULL,
     r = r,
     lambda.a = lambda.a,
     S.a = S.a,
     lambda.b = lambda.b,
     S.b = S.b,
     warm.start = warm,
     trace.it = F,
     thresh = thresh,
     final.svd = T,
     maxit = maxit
    )
   
   #--------------------------------------------------------------
   # predicting validation set and xbetas for next fit:
   XbetaValid = suvC(X%*%fiti$ub, fiti$vb %*% fiti$db, virow, vpcol)
   MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), virow, vpcol)
   #--------------------------------------------
   err = error_function(MValid + XbetaValid, y_valid)
   #rank <- sum(round(fiti$d, 4) > 0)
   # newly added, to be removed later
   var_explained = fiti$d ^ 2 / sum(fiti$d ^ 2)
   cum_var = cumsum(var_explained)
   rank  <- which(cum_var >= pct)[1]
   
   warm <- fiti # warm start for next
   #---------------------------------------------------------------------
   if (trace == TRUE)
    print(sprintf(
     paste0(
      "%2d lambda.M=%9.5g, rank.max = %d  ==>",
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
    best_fit$lambda.M = lamseq[i]
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
    CASMC2_fit2(
     y = y,
     X = X,
     J = best_fit$rank.max,
     lambda.M = best_fit$lambda.M,
     r = r,
     svdH = NULL,
     lambda.beta = lambda.beta,
     lambda.a = lambda.a,
     S.a = S.a,
     lambda.b = lambda.b,
     S.b = S.b,
     warm.start = best_fit$fit,
     trace.it = F,
     thresh = thresh,
     final.svd = T,
     maxit = maxit
    )
   
  }
  best_fit$lambda.beta = lambda.beta
  best_fit$lambda.a = lambda.a
  best_fit$lambda.b = lambda.b
  return(best_fit)
 }
#---