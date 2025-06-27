#' Soft‐thresholding operator
#'
#' Apply the element‐wise soft‐thresholding \(S_\lambda(x)=\mathrm{sign}(x)\max(|x|-\lambda,0)\)
#'
#' @param B Numeric matrix
#' @param lambda Non‐negative scalar threshold
#' @return Numeric matrix of same dimensions as `B`
#' @export
soft_threshold <- function(B, lambda) {
  sign(B) * pmax(abs(B) - lambda, 0)
}

#' Covariate‐Adjusted‐Sparse‐Matrix‐Completion Fit
#'
#' Fit low‐rank plus covariate model
#'
#' @inheritParams soft_threshold
#' @param y CsparseMatrix of class “Incomplete”
#' @param X Covariate matrix
#' @param J Maximum rank of low‐rank component (default = 2)
#' @param lambda_M Nuclear‐norm penalty for low‐rank component (default = 0)
#' @param lambda_beta Lasso penalty for β (default = 0)
#' @param S_a,lambda_a Optional similarity matrix + weight for rows (A)
#' @param S_b,lambda_b Optional similarity matrix + weight for cols (B)
#' @param normalized_laplacian Logical, whether to normalize Laplacian (default = TRUE)
#' @param maxit Maximum number of iterations (default = 300)
#' @param thresh Convergence threshold on Frobenius‐change (default = 1e-5)
#' @param trace Logical, if TRUE prints progress (default = FALSE)
#' @param warm_start Optional warm‐start list with components `u`,`d`,`v`,`beta`
#' @param final_svd Logical, whether to do final SVD refinement (default = TRUE)
#' @param min_eigv Minimum eigenvalue floor (default = 0)
#' @return A list with  
#'   - `u,d,v`: factors of low‐rank component  
#'   - `beta`: covariate coefficients  
#'   - hyperparameters and `n_iter`  
#' @export
CAMC3_fit <- CAMC_Lasso_fit <- function(
    y,
    X,
    J                 = 2,
    lambda_M          = 1,
    lambda_beta       = 0,
    S_a               = NULL,
    lambda_a          = 0,
    S_b               = NULL,
    lambda_b          = 0,
    normalized_laplacian = TRUE,
    maxit             = 300,
    thresh            = 1e-5,
    trace             = FALSE,
    warm_start        = NULL,
    final_svd         = TRUE,
    min_eigv          = 0
) {
  ## 1. Input checks & setup
  stopifnot(inherits(y, "CsparseMatrix"))
  i_row  <- y@i
  p_col  <- y@p
  y_obs  <- y@x
  dims   <- dim(y)
  n_rows <- dims[1]
  n_cols <- dims[2]
  if (trace) nz <- nnzero(y, na.counted = TRUE)
  
  ## 2. Laplacian reg. if requested
  laplace_a <- FALSE
  laplace_b <- FALSE
  if (!is.null(S_a) && lambda_a > 0) {
    laplace_a <- TRUE
    L_a <- utils$computeLaplacian(S_a,
                                  normalized = normalized_laplacian) * lambda_a
  }
  if (!is.null(S_b) && lambda_b > 0) {
    laplace_b <- TRUE
    L_b <- utils$computeLaplacian(S_b,
                                  normalized = normalized_laplacian) * lambda_b
  }
  
  ## 3. Warm‐start or initialize
  clean.warm.start(warm_start)
  if (!is.null(warm_start)) {
    if (!all(match(c("u","d","v","beta"),
                   names(warm_start), 0) > 0)) {
      stop("warm_start missing components: u, d, v, beta")
    }
    beta <- warm_start$beta
    ws   <- utils$prepare.M.warm.start(
      warm_start, J, n_rows, n_cols, min_eigv)
    U    <- ws$U;    V    <- ws$V;    Dsq  <- ws$Dsq
  } else {
    beta <- as.matrix(crossprod(X, y))
    y@x  <- y_obs - suvC(X, t(beta), i_row, p_col)
    init <- utils$svdopt(naive_MC(as.matrix(y)),
                         J, n_rows, n_cols, FALSE, FALSE)
    U    <- init$u;  Dsq  <- init$d;  V    <- init$v
  }
  
  ## 4. Main loop
  ratio <- Inf
  iter  <- 0
  while (ratio > thresh && iter < maxit) {
    iter <- iter + 1
    U_old  <- U;  V_old  <- V;  D_old <- Dsq
    
    # 4.1 Update residuals
    VDsq    <- t(Dsq * t(V))
    M_obs   <- suvC(U, VDsq, i_row, p_col)
    if (iter == 1) {
      xb_obs <- suvC(X, t(beta), i_row, p_col)
      y@x   <- y_obs - M_obs - xb_obs
    }
    
    # 4.2 Update beta via soft‐threshold
    beta   <- soft_threshold(
      as.matrix(crossprod(X, y) + beta),
      lambda_beta
    )
    xb_obs <- suvC(X, t(beta), i_row, p_col)
    y@x   <- y_obs - M_obs - xb_obs
    
    # 4.3 Update B → (V, Dsq, U)
    B_mat  <- crossprod(y, U) + VDsq
    if (laplace_b) B_mat <- B_mat - t(V) %*% L_b
    D_star <- Dsq / (Dsq + lambda_M)
    B_s    <- as.matrix(B_mat %*% diag(D_star))
    Bs     <- utils$svd_small_nc(B_s, FALSE, p = J)
    V      <- Bs$u;  Dsq <- Bs$d; U <- U %*% Bs$v
    
    # 4.4 Update A → (U, Dsq, V)
    A_mat  <- y %*% V + t(Dsq * t(U))
    if (laplace_a) A_mat <- A_mat - L_a %*% U
    D_star <- Dsq / (Dsq + lambda_M)
    A_s    <- as.matrix(A_mat %*% diag(D_star))
    As     <- utils$svd_small_nc(A_s, FALSE, p = J)
    U      <- As$u;  Dsq <- As$d; V <- V %*% As$v
    
    # 4.5 Convergence check
    ratio <- utils$Frob(U_old, D_old, V_old, U, Dsq, V)
    if (trace) {
      obj <- (0.5 * sum(y@x^2) + lambda_M * sum(Dsq)) / nz
      cat(iter, " obj=", round(obj,5), " ratio=", ratio, "\n")
    }
  }
  if (iter == maxit && trace)
    warning("Did not converge in ", maxit, " iterations.")
  
  ## 5. Final SVD refinement
  if (final_svd) {
    A_mat <- y %*% V + t(Dsq * t(U))
    if (laplace_a) A_mat <- A_mat - L_a %*% U
    A_fin <- t(t(A_mat) * (Dsq / (Dsq + lambda_M)))
    fin   <- utils$svd_small_nc(as.matrix(A_fin), FALSE, p = J)
    U     <- fin$u; V <- V %*% fin$v
    Dsq   <- pmax(fin$d - lambda_M, 0)
    if (trace) {
      M_obs <- suvC(t(Dsq * t(U)), V, i_row, p_col)
      y@x  <- y_obs - M_obs - xb_obs
      obj  <- (0.5 * sum(y@x^2) + lambda_M * sum(Dsq)) / nz
      cat("final SVD obj=", round(obj,5), "\n",
          "iterations=", iter, "\n")
    }
  }
  
  ## 6. Trim effective rank and return
  J_eff <- min(max(1, sum(Dsq > min_eigv)), J)
  list(
    u           = U[, seq_len(J_eff), drop = FALSE],
    d           = Dsq[seq_len(J_eff)],
    v           = V[, seq_len(J_eff), drop = FALSE],
    beta        = beta,
    lambda_M    = lambda_M,
    lambda_beta = lambda_beta,
    J           = J_eff,
    lambda_a    = lambda_a,
    lambda_b    = lambda_b,
    n_iter      = iter
  )
}
