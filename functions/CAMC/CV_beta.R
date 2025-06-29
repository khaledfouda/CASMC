#' Covariate‐Adjusted‐Sparse‐Matrix‐Completion Lasso Cross‐Validation
#'
#' Holdout‐set CV over a grid of \eqn{λ_β}, calling \code{CAMC3_cv_M} under the hood.
#'
#' @param y_train CsparseMatrix (Incomplete) of training responses.
#' @param X Covariate matrix.
#' @param y_valid Numeric vector of held‐out responses.
#' @param W_valid Indicator matrix (0 = validation, 1 = training).
#' @param y Optional full CsparseMatrix (train + valid) for final refit.
#' @param hpar List of hyperparameters (see \code{CAMC_Lasso_hparams}).
#' @param error_function Function to compute validation error (default = RMSE).
#' @param thresh Convergence tolerance (default = 1e-6).
#' @param maxit Maximum iterations per fit (default = 100).
#' @param verbose Integer: 0 = silent, 1 = per‐fit summary, 2 = verbose (passed to inner fits).
#' @param max_cores Maximum parallel workers (default = 8).
#' @param seed Optional RNG seed.
#'
#' @return A list containing the best fit object, its error, hyperparameters, and init params.
#' @export
CAMC_Lasso_cv <- CAMC_cv <- function(
    y_train,
    X,
    y_valid,
    W_valid,
    y               = NULL,
    hpar            = CAMC_Lasso_hparams,
    error_function  = utils$error_metric$rmse,
    thresh          = 1e-6,
    maxit           = 100,
    verbose         = 0,
    max_cores       = 8,
    seed            = NULL
) {
  ##-------------------------------------------------------------------------------
  ## 1. Reproducibility & basic checks
  if (!is.null(seed)) set.seed(seed)
  stopifnot(inherits(y_train, "CsparseMatrix"))
  
  ##-------------------------------------------------------------------------------
  ## 2. Get an upperbound to lambda
  if (is.null(hpar$beta$lambda_max) | !(is.numeric(hpar$beta$lambda_max)))
    hpar$beta$lambda_max <- find_lasso_max_param(
      y_train = y_train,
      X       = X,
      y_valid = y_valid,
      W_valid = W_valid,
      y       = y,
      maxit   = 100,
      verbose = verbose,
    )
    
  ## 3. Build λ_β grid
  lambda_beta_grid <- seq(
    from       = hpar$beta$lambda_max,
    to         = .Machine$double.eps,
    length.out = hpar$beta$n.lambda
  )
  
  ##-------------------------------------------------------------------------------
  ## 4. Parallel setup
  num_cores <- min(detectCores(logical = FALSE) - 1,
                   length(lambda_beta_grid))
  if (num_cores > max_cores) 
    num_cores <- min(max_cores, ceiling(num_cores / 2))
  
  if(verbose > 0)
    message("Running on ", num_cores, " cores.")
  
  ##-------------------------------------------------------------------------------
  ## 5. CV loop: fit CAMC3_cv_M for each λ_β in parallel
  inner_trace <- (verbose >= 2)
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  plan(multisession, workers = num_cores)
  results <- future_lapply(
    lambda_beta_grid,
    function(lambda_beta) {
      tryCatch(
        CAMC3_cv_M(
          y_train        = y_train,
          X              = X,
          y_valid        = y_valid,
          W_valid        = W_valid,
          y              = y,
          hpar           = hpar,
          lambda_beta    = lambda_beta,
          error_function = error_function,
          thresh         = thresh,
          maxit          = maxit,
          trace          = inner_trace,
          seed           = seed
        ),
        error = function(e) {
          list(
            error_message = e,
            error         = Inf,
            lambda_beta   = lambda_beta
          )
        }
      )
    }, future.seed=TRUE, future.globals = TRUE, future.packages = packages
  )
  ##-------------------------------------------------------------------------------
  ## 6. Report any fit errors
  for (res in results) {
    if (!is.null(res$error_message)) {
      message(
        "Error at λ_β=", round(res$lambda_beta, 3), 
        ": ", res$error_message
      )
    }
  }
  
  ##-------------------------------------------------------------------------------
  ## 7. Select the best fit
  errors    <- vapply(results, `[[`, numeric(1), "error")
  best_idx  <- which.min(errors)
  best_fit  <- results[[best_idx]]
  
  ##-------------------------------------------------------------------------------
  ## 8. Optional printing
  if (verbose >= 1) {
    for (res in results) {
      message(sprintf(
        "<< λ_β=%.4g | sparsity=%.2f | err=%.5f | iters=%d | rank_M=%d | λ_M=%.4g >>",
        res$lambda_beta,
        sum(res$fit$beta == 0) / length(res$fit$beta),
        res$error,
        res$fit$n_iter,
        res$fit$J,
        res$fit$lambda_M
      ))
    }
    message(sprintf(
      "<< Best fit >> λ_β=%.4g | sparsity=%.2f | err=%.5f | iters=%d | rank_M=%d | λ_M=%.4g >>",
      best_fit$lambda_beta,
      sum(best_fit$fit$beta == 0) / length(best_fit$fit$beta),
      best_fit$error,
      best_fit$fit$n_iter,
      best_fit$fit$J,
      best_fit$fit$lambda_M
    ))
  }
  
  ##-------------------------------------------------------------------------------
  ## 9. Attach metadata & return
  best_fit$init_hpar <- hpar
  
  best_fit
}