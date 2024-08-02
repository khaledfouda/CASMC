#' Covariate-Adjusted-Sparse-Matrix-completion
#' Fit function
#'
#' @param y A sparse matrix of class Incomplete.
#' @param X covariate matrix
#' @param svdH (optional) A list consisting of the SVD of the hat matrix. see reduced_hat_decomp function.
#' @param Xterms (optional) A list of terms computed using GetXterms function.
#' @param J (hyperparameter). The maximum rank of the low-rank matrix M = AB. Default is 2
#' @param lambda (hyperparameter) The L2 regularization parameter for A and B
#' @param r (hyperparameter, optional) The rank of covariate effects (beta). The rank will be the same
#'                                      as the covariate matrix if r not provided
#' @return A list of u,d,v of M, Beta, and a vector of the observed Xbeta
#' @examples
#'  CASMC_fit(y,X,J=5)
#' @export
#'
CASMC1_fit <-
  CASMC_Ridge_fit <-
  function(y,
           X,
           #svdH = NULL,
           Xterms = NULL,
           J = 2,
           #r = NULL,
           lambda = 0,
           lambda.beta = 0, # not needed if Xterms is provided
           # similarity matrix for A
           S.a = NULL,
           lambda.a = 0,
           # similarity matrix for B
           S.b = NULL,
           lambda.b = 0,
           maxit = 100,
           thresh = 1e-05,
           trace.it = FALSE,
           warm.start = NULL,
           final.svd = TRUE,
           min_eigv = 1e-4) {
    stopifnot(inherits(y, "dgCMatrix"))
    irow = y@i
    pcol = y@p
    n <- dim(y)
    m <- n[2]
    n <- n[1]
    k <- ncol(X)
    if (trace.it)
      nz = nnzero(y, na.counted = TRUE)
    #-------------------------------
    laplace.a <- laplace.b <- F
    if (!is.null(S.a) & lambda.a > 0) {
      laplace.a = T
      L.a = utils$computeLaplacian(S.a, normalized = F) * lambda.a
    }
    if (!is.null(S.b) & lambda.b > 0) {
      laplace.b = T
      L.b = utils$computeLaplacian(S.b, normalized = F) * lambda.b
    }
    #--------------------------------
    # if svdH is not given but X is given. only needed if warm.start is not provided
    if (is.null(warm.start))
      svdH = utils$reduced_hat_decomp.H(X)
    #---------------------------------------------------
    # if Xterms are not provided but X is given.
    if (is.null(Xterms))
      Xterms = utils$GetXterms(X, lambda.beta)
    #---------------------------------------------------
    # warm start or initialize (naive or random)
    warm = FALSE
    clean.warm.start(warm.start)
    if (!is.null(warm.start)) {
      #must have u,d and v components
      if (!all(match(c("u", "d", "v", "beta"), names(warm.start), 0) >
               0))
        stop("warm.start does not have components u, d and v")
      warm = TRUE
      beta = warm.start$beta
      warm.out <- utils$prepare.M.warm.start(warm.start, J, n, m, min_eigv)
      U = warm.out$U
      V = warm.out$V
      Dsq = warm.out$Dsq
    } else{
      # initialize. Warm start is not provided
      
      Y_naive = naive_MC(as.matrix(y))
      beta = Xterms$X1 %*% Y_naive
      Xbeta <- X %*% beta   #svdH$u %*% (svdH$v  %*% Y_naive)
      M <- Y_naive - Xbeta
      M <- utils$svdopt(as.matrix(M), J, n, m, F, F)
      U = M$u
      V = M$v
      Dsq = M$d
      #----------------------
      # initialization for beta = X^-1 Y
      # comment for later: shouldn't be X^-1 H Y??
      #beta = as.matrix(ginv(X) %*% Xbeta)
      Y_naive <- Xbeta <- M <- NULL
      #---------------------------------------------------------------
      
    }
    #----------------------------------------
    yobs <- y@x
    ratio <- 1
    iter <- 0
    #----------------------------------------
    while ((ratio > thresh) & (iter < maxit)) {
      iter <- iter + 1
      U.old = U
      V.old = V
      Dsq.old = Dsq
      #----------------------------------------------
      # Part 1: Update Beta while A and B are fixed
      # updates xbeta.obs and Beta.
      VDsq = t(Dsq * t(V))
      M_obs = suvC(U, VDsq, irow, pcol)
      if (iter == 1)
        y@x = yobs - M_obs #- suvC(X, t(beta), irow, pcol)
      beta = as.matrix(Xterms$X1 %*% y + (Xterms$X2 %*%  beta))
      xbeta.obs <- suvC(X, t(beta), irow, pcol)
      y@x = yobs - M_obs - xbeta.obs
      ##--------------------------------------------
      # part 2: Update B while A and beta are fixed
      # updates U, Dsq, V
      B = t(U) %*% y + t(VDsq)
      if (laplace.b)
        B = B - t(V)  %*% L.b
      B = t((B) * (Dsq / (Dsq + lambda)))
      Bsvd = utils$svd_small_nc(as.matrix(B), FALSE, p = J)  
      V = Bsvd$u
      Dsq = Bsvd$d
      U = U %*% (Bsvd$v)
      #-------------------------------------------------------------
      # part 3: Update A while B and beta are fixed
      # updates U, Dsq, V
      A = (y %*% V) + t(Dsq * t(U))
      if (laplace.a)
        A = A - L.a %*% U
      A = t(t(A) * (Dsq / (Dsq + lambda)))
      Asvd =  utils$svd_small_nc(as.matrix(A), FALSE, p = J)
      U = Asvd$u
      Dsq = Asvd$d
      V = V %*% (Asvd$v)
      #------------------------------------------------------------------------------
      ratio =  utils$Frob(U.old, Dsq.old, V.old, U, Dsq, V)
      #------------------------------------------------------------------------------
      if (trace.it) {
        obj = (.5 * sum(y@x ^ 2) + lambda * sum(Dsq)) / nz
        cat(iter, ":", "obj", format(round(obj, 5)), "ratio", ratio, "\n")
      }
      #-----------------------------------------------------------------------------------
      
    }
    if (iter == maxit &
        trace.it)
      warning(
        paste(
          "Convergence not achieved by",
          maxit,
          "iterations. Consider increasing the number of iterations."
        )
      )
    
    # one final fit for one of the parameters (A) has proved to improve the performance significantly.
    if (final.svd) {
      #---- update A
      A = (y %*% V) + t(Dsq * t(U))
      if (laplace.a)
        A = A - L.a %*% U
      A = t(t(A) * (Dsq / (Dsq + lambda)))
      Asvd =  utils$svd_small_nc(as.matrix(A), FALSE, p = J)
      U = Asvd$u
      Dsq = Asvd$d
      V = V %*% (Asvd$v)
      Dsq = pmax(Dsq - lambda, 0)
      #------------------
      if (trace.it) {
        M_obs = suvC(t(Dsq * t(U)), V, irow, pcol)
        y@x = yobs - M_obs - xbeta.obs
        obj = (.5 * sum(y@x ^ 2) + lambda * sum(Dsq)) / nz
        cat("final SVD:", "obj", format(round(obj, 5)), "\n")
        cat("Number of Iterations for covergence is ", iter, "\n")
      }
    }
    #-------------------------------------------------------
    # trim in case we reduce the rank of M to be smaller than J.
    J = min(max(1,sum(Dsq > min_eigv)), J)
    out = list(
      u = U[, seq(J), drop = FALSE],
      d = Dsq[seq(J)],
      v = V[, seq(J), drop = FALSE],
      lambda = lambda,
      J = J,
      n_iter = iter,
      beta = beta
    )
    out
  }
