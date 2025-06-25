#' Generate simulated data for CAMC
#'
#' Simulate covariates X, low-rank structure M orthogonal to X, coefficients β, 
#' and observed matrix Y = (Xβ + M + noise) masked by missingness.
#'
#' @param n Number of samples (rows) (default = 300)
#' @param m Number of features (columns) (default = 400)
#' @param r Rank of low-rank component M (default = 10)
#' @param k Number of covariates (default = 10)
#' @param missing_prob Probability an entry is missing (default = 0.8)
#' @param collinear Logical; if TRUE, make columns 1 and 2 of X nearly collinear (default = FALSE)
#' @param half_discrete Logical; if TRUE, randomly binarize half of X’s columns (default = FALSE)
#' @param informative_cov_prop Proportion of β entries that are nonzero (default = 1)
#' @param prepare_for_fitting Logical; if TRUE, also return a train/test split for fitting (default = FALSE)
#' @param mar_sparse Logical; if TRUE and informative_cov_prop<1, zero a random subset of β entries (default = FALSE)
#' @param mv_beta Logical; if TRUE, draw β rows from a multivariate normal; otherwise uniform (default = TRUE)
#' @param seed Optional random seed for reproducibility
#' @return A list with components:
#'   - `O`: noise-free matrix Xβ + M  
#'   - `W`: missingness mask (1 = observed, 0 = missing)  
#'   - `X`: covariate matrix  
#'   - `Y`: observed matrix, (O + noise) * W  
#'   - `beta`: coefficient matrix (k × m)  
#'   - `M`: low-rank component (n × m)  
#'   - `rank`: numeric rank of O  
#'   - `fit_data` (optional): list with `train`, `valid`, `Y_full`, `W_valid` if `prepare_for_fitting = TRUE`  
#' @export
generate_simulated_data <- function(
    n                     = 300,
    m                     = 400,
    r                     = 10,
    k                     = 10,
    missing_prob          = 0.8,
    collinear             = FALSE,
    half_discrete         = FALSE,
    informative_cov_prop  = 1,
    prepare_for_fitting   = FALSE,
    mar_sparse            = FALSE,
    mv_beta               = TRUE,
    seed                  = NULL
) {
  ## 1. Reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ## 2. Simulate covariates X
  X <- matrix(rnorm(n * k), nrow = n, ncol = k)
  if (half_discrete) {
    n_disc     <- ceiling(k / 2)
    disc_cols  <- sample(seq_len(k), n_disc)
    X[, disc_cols] <- rbinom(n * n_disc, size = 1, prob = 0.3)
  }
  
  if (collinear) {
    # make column 2 nearly a copy of column 1
    X[, 2] <- X[, 1] + rnorm(n, mean = 0, sd = 1.5)
  }
  
  ## 3. Simulate coefficient matrix beta (k × m)
  if (mv_beta) {
    beta_means <- runif(k, 1, 3) * sample(c(-1, 1), k, replace = TRUE)
    beta_vars  <- runif(k, 0.5, 1)^2
    beta <- t(MASS::mvrnorm(
      n   = m,
      mu  = beta_means,
      Sigma = diag(beta_vars, nrow = k)
    ))
  } else {
    beta <- matrix(
      runif(k * m, 1, 2),
      nrow = k,
      ncol = m
    )
  }
  
  ## 4. Simulate low-rank structure M orthogonal to X
  U      <- matrix(runif(n * r), nrow = n,       ncol = r)
  V      <- matrix(runif(r * m), nrow = r,       ncol = m)
  P_X    <- X %*% solve(crossprod(X)) %*% t(X)       # projection onto col(X)
  P_bar  <- diag(n) - P_X                            # projection orthogonal to col(X)
  M      <- P_bar %*% U %*% V
  
  ## 5. Generate missingness mask W (1 = observed, 0 = missing)
  W <- matrix(
    rbinom(n * m, size = 1, prob = 1 - missing_prob),
    nrow = n,
    ncol = m
  )
  
  ## 6. Enforce proportion of informative covariates
  if (informative_cov_prop < 1) {
    total_beta <- length(beta)
    if (mar_sparse) {
      to_zero <- sample(
        seq_len(total_beta),
        size = round((1 - informative_cov_prop) * total_beta)
      )
      beta[to_zero] <- 0
    } else {
      n_keep <- round(informative_cov_prop * k)
      if (n_keep <= 0) {
        beta[] <- 0
      } else if (n_keep < k) {
        remove_idx <- sample(seq_len(k), k - n_keep)
        beta[remove_idx, ] <- 0
      }
    }
  }
  
  ## 7. Construct observed data Y with additive Gaussian noise
  O        <- X %*% beta + M
  noise_sd <- 1
  E        <- matrix(
    rnorm(n * m, mean = 0, sd = noise_sd),
    nrow = n,
    ncol = m
  )
  Y        <- (O + E) * W
  rank_O   <- qr(O)$rank
  
  ## 8. Optional train/test split for model fitting
  fit_data <- NULL
  if (prepare_for_fitting) {
    W_valid      <- utils$MC_train_test_split(W, testp = 0.2)
    Y_train      <- Y * W_valid
    Xq <- qr(X)
    fit_data     <- list(
      train   = utils$to_incomplete(Y_train),
      valid   = Y[W_valid == 0],
      Y_full  = utils$to_incomplete(Y),
      W_valid = W_valid,
      Xq      = qr.Q(Xq),
      Xr      = qr.R(Xq),
      Rbeta   = qr.R(Xq) %*% beta
    )
  }
  
  ## 9. Return all components
  list(
    O         = O,
    W         = W,
    X         = X,
    Y         = Y,
    beta      = beta,
    M         = M,
    rank      = rank_O,
    fit_data  = fit_data
  )
}