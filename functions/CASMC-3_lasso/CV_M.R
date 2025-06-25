
#---------------------------------------------------------------------------------------------
#' Covariate‐Adjusted‐Sparse‐Matrix‐Completion Cross‐Validation
#'
#' Holdout cross‐validation for CAMC3, fitting on `y_train` and evaluating on `y_valid`.
#'
#' @param y_train A dgCMatrix (Incomplete) training response.
#' @param X Covariate matrix.
#' @param y_valid Numeric vector of held‐out responses.
#' @param W_valid Indicator matrix (0 = validation, 1 = train); used to extract validation indices.
#' @param y Optional full dgCMatrix (train + validation).  If provided, a final fit on the full data is done.
#' @param lambda_beta L2‐penalty for β (default = machine epsilon).
#' @param hpar List of hyperparameters (must contain `$M`, `$beta`, and `$laplacian` sub‐lists).
#' @param error_function Function to measure validation error (default = `utils$error_metric$rmse`).
#' @param thresh Convergence threshold (default = 1e‐6).
#' @param maxit Maximum number of iterations per fit (default = 100).
#' @param trace If `TRUE`, prints progress messages.
#' @param warm Optional warm‐start list (with components `u, d, v, beta`).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A list with elements:
#'   - `error`: best validation error  
#'   - `rank_M`: rank of M at best lambda  
#'   - `lambda_M`: best nuclear‐penalty lambda  
#'   - `rank_max`: rank.max used in best fit  
#'   - `fit`: the corresponding `CAMC3_fit` object  
#'   - `iter`: index in the lambda‐sequence  
#'   - plus the final `lambda_beta`, `lambda_a`, `lambda_b`, `learning.rate`  
#' @export
CAMC3_cv_M <- function(
    y_train,
    X,
    y_valid,
    W_valid,
    y            = NULL,
    lambda_beta  = .Machine$double.eps,
    hpar,
    error_function = utils$error_metric$rmse,
    thresh         = 1e-6,
    maxit          = 100,
    trace          = FALSE,
    warm           = NULL,
    seed           = NULL
) {
  ##------------------------------------------------------------------------------
  ## 1. Reproducibility & Input checks
  if (!is.null(seed)) set.seed(seed)
  stopifnot(inherits(y_train, "dgCMatrix"))
  
  ##------------------------------------------------------------------------------
  ## 2. Build lambda‐sequence for nuclear penalty
  if (is.null(hpar$M$lambda_init)) {
    hpar$M$lambda_init <- utils$lambdaM.max(
      y_train,
      utils$reduced_hat_decomp.H(X)
    ) * hpar$M$lambda_factor
  }
  lambda_seq <- seq(
    from      = hpar$M$lambda_init,
    to        = .Machine$double.eps,
    length.out= hpar$M$n.lambda
  )
  
  ##------------------------------------------------------------------------------
  ## 3. Extract validation indices from W_valid
  W_valid[W_valid == 1] <- NA
  W_valid[W_valid == 0] <- 1
  W_valid <- as(W_valid, "Incomplete")
  vi_row   <- W_valid@i
  vp_col   <- W_valid@p
  rm(W_valid)
  
  ##------------------------------------------------------------------------------
  ## 4. Initialize CV loop
  rank_max          <- hpar$M$rank.init
  best_fit          <- list(error = Inf)
  no_improve_count  <- 0
  
  for (i in seq_along(lambda_seq)) {
    ## 4.1 Fit on training data
    fit_i <- CAMC3_fit(
      y             = y_train,
      X             = X,
      J             = rank_max,
      lambda_M      = lambda_seq[i],
      lambda_beta   = lambda_beta,
      lambda_a      = hpar$laplacian$lambda_a,
      S_a           = hpar$laplacian$S.a,
      lambda_b      = hpar$laplacian$lambda_b,
      S_b           = hpar$laplacian$S.b,
      warm_start    = warm,
      trace         = FALSE,
      thresh        = thresh,
      maxit         = maxit
    )
    
    ## 4.2 Predict & compute validation error
    xbeta_valid <- suvC(X, t(fit_i$beta), vi_row, vp_col)
    m_valid     <- suvC(
      fit_i$u,
      t(fit_i$d * t(fit_i$v)),
      vi_row,
      vp_col
    )
    error_val   <- error_function(m_valid + xbeta_valid, y_valid)
    current_rank<- sum(round(fit_i$d, 4) > 0)
    
    # prepare warm start for next lambda
    warm <- fit_i
    
    if (trace) {
      message(sprintf(
        "%2d lambda=%.4g | rank_max=%d => rank=%d | err=%.5f | iters=%d",
        i,
        lambda_seq[i],
        rank_max,
        current_rank,
        error_val,
        fit_i$n_iter
      ))
    }
    
    ## 4.3 Track best model & early stopping
    if (error_val < best_fit$error) {
      best_fit <- list(
        error     = error_val,
        rank_M    = current_rank,
        lambda_M  = lambda_seq[i],
        rank_max  = rank_max,
        fit       = fit_i,
        iter      = i
      )
      no_improve_count <- 0
    } else {
      no_improve_count <- no_improve_count + 1
    }
    
    if (no_improve_count >= hpar$M$early.stopping) {
      if (trace) {
        message(
          sprintf(
            "Early stopping: no improvement in last %d lambda’s.",
            no_improve_count
          )
        )
      }
      break
    }
    
    # update rank_max for next iteration
    rank_max <- min(
      current_rank + hpar$M$rank.step,
      hpar$M$rank.limit
    )
  }
  
  ##------------------------------------------------------------------------------
  ## 5. (Optional) Final fit on full data if `y` is provided
  if (!is.null(y)) {
    stopifnot(inherits(y, "dgCMatrix"))
    best_fit$fit <- CAMC3_fit(
      y             = y,
      X             = X,
      J             = best_fit$rank_max,
      lambda_M      = best_fit$lambda_M,
      lambda_beta   = lambda_beta,
      lambda_a      = hpar$laplacian$lambda_a,
      S_a           = hpar$laplacian$S.a,
      lambda_b      = hpar$laplacian$lambda_b,
      S_b           = hpar$laplacian$S.b,
      warm_start    = best_fit$fit,
      trace         = FALSE,
      thresh        = thresh,
      maxit         = maxit
    )
  }
  
  ##------------------------------------------------------------------------------
  ## 6. Record final hyperparameters & return
  best_fit$lambda_beta   <- lambda_beta
  best_fit$lambda_a      <- hpar$laplacian$lambda_a
  best_fit$lambda_b      <- hpar$laplacian$lambda_b
  best_fit$learning.rate <- hpar$beta$learning.rate
  
  best_fit
}
