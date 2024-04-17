# this one incorporates SZIRCI into the model:
# Option 2: get rid of the Xbeta and use the new estimates!
GetXterms <- function(X) {
  svdX = fast.svd(X)
  Ux = svdX$u
  Vx = svdX$d * t(svdX$v)
  X0 = ginv(t(Vx) %*% Vx) %*% t(Vx)
  X1 = X0 %*% t(Ux)
  X2 = X0 %*% Vx
  list(X1 = X1, X2 = X2)
}

#' Sum Two Numbers
#'
#' This function takes two numeric inputs and returns their sum.
#' @param y A sparse matrix of class Incomplete.
#' @param svdH A list consisting of the SVD of the hat matrix. see reduced_hat_decomp function.
#' @param X optional
#' @param J Hyperparameter. The maximum rank of the low-rank matrix M = AB
#' @param Xterms ...
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(1, 3)
#' @export
#'
CASMC_fit <-
  function(y,
           X = NULL,
           svdH = NULL,
           J = 2,
           thresh = 1e-05,
           lambda = 0,
           maxit = 100,
           trace.it = FALSE,
           warm.start = NULL,
           Xterms = NULL,
           final.svd = TRUE,
           init = "naive") {
    stopifnot(inherits(y, "dgCMatrix"))
    irow = y@i
    pcol = y@p
    n <- dim(y)
    m <- n[2]
    n <- n[1]
    nz = nnzero(y, na.counted = TRUE)
    initialize_beta = FALSE
    beta.obs = NULL
    #-------------------------------
    # if svdH is not given but X is given.
    if (is.null(svdH)) {
      stopifnot(!is.null(X))
      svdH = reduced_hat_decomp(X, 1e-2)
      J_H = svdH$rank
      svdH = svdH$svdH
      if (trace.it)
        print(paste("Rank of H is ", J_H))
    }
    #---------------------------------------------------
    # if Xterms are not provided but X is given.
    if (is.null(Xterms)) {
      stopifnot(!is.null(X))
      Xterms = GetXterms(X)
    }
    X1 = Xterms$X1
    X2 = Xterms$X2
    #---------------------------------------------------
    # warm start or initialize (naive or random)
    warm = FALSE
    clean.warm.start(warm.start)
    if (!is.null(warm.start)) {
      #must have u,d and v components
      if (!all(match(c("u", "d", "v", "xbeta.obs", "beta"), names(warm.start), 0) >
               0))
        stop("warm.start does not have components u, d and v")
      warm = TRUE
      D = warm.start$d
      JD = sum(D > 0)
      #J = JD
      beta = as.matrix(warm.start$beta)
      Beta = fast.svd(beta)
      #xbeta.obs = warm.start$xbeta.obs
      xbeta.obs <-
        suvC(X %*% Beta$v, t(Beta$d * t(Beta$u)), irow, pcol)
      if (JD >= J) {
        U = warm.start$u[, seq(J), drop = FALSE]
        V = warm.start$v[, seq(J), drop = FALSE]
        Dsq = D[seq(J)]
      } else{
        Dsq = c(D, rep(D[JD], J - JD))
        Ja = J - JD
        U = warm.start$u
        Ua = matrix(rnorm(n * Ja), n, Ja)
        Ua = Ua - U %*% (t(U) %*% Ua)
        Ua = fast.svd(Ua)$u
        U = cbind(U, Ua)
        V = cbind(warm.start$v, matrix(0, m, Ja))
      }
    } else{
      # initialize. Warm start is not provided
      stopifnot(init %in% c("random", "naive"))
      if (init == "random") {
        V = matrix(0, m, J)
        U = matrix(rnorm(n * J), n, J)
        U = fast.svd(U)$u
        Dsq = rep(1, J)# we call it Dsq because A=UD and B=VDsq and AB'=U Dsq V^T
        xbeta.obs = suvC(svdH$u, t(as.matrix(svdH$v %*% y)), irow, pcol)
        stop("Beta is random initialization isn't implemented yet.")
      } else if (init == "naive") {
        Y_naive = as.matrix(y)
        yobs = Y_naive != 0
        Y_naive = naive_MC(Y_naive)
        naive_fit <-  svdH$u %*% (svdH$v  %*% Y_naive)
        # Xbeta = H Y
        xbeta.obs <- naive_fit[yobs]
        # M = (I-H) Y
        naive_fit <- Y_naive - naive_fit
        naive_fit <- propack.svd(as.matrix(naive_fit), J)
        U = naive_fit$u
        V = naive_fit$v
        Dsq = naive_fit$d
        #----------------------
        # initialization for SZIRCI
        beta = t(ginv(X) %*% Y_naive) # B = (X^-1 Y)'
        Beta = fast.svd(beta)
        # xbeta.obs <- suvC(X %*% Beta$v, t(Beta$d * t(Beta$u)), irow, pcol)
        #---------------------------------------------------------------
      }
    }
    
    #----------------------------------------
    S <- y
    ratio <- 1
    iter <- 0
    counter = 0
    best_score = Inf
    best_iter = NA
    #-----------------------------------------------------------
    
    #--------------------------------------------------------------
    #####
    # update Beta once at the beginning if parameters were initialized
    if (is.null(warm.start)) {
      VDsq = t(Dsq * t(V))
      HU = svdH$u %*% (svdH$v %*% U)
      xbeta.obs = xbeta.obs +
        suvC(svdH$u, t(as.matrix(svdH$v %*% S)), irow, pcol) +
        suvC(HU, VDsq, irow, pcol)
    }
    # initial M and S (sparse Y - M - Xbeta )
    M_obs = suvC(U, VDsq, irow, pcol)
    S@x = y@x - M_obs - xbeta.obs
    #----------------------------------------
    while ((ratio > thresh) & (iter < maxit)) {
      iter <- iter + 1
      U.old = U
      V.old = V
      Dsq.old = Dsq
      #----------------------------------------------
      # Part 1: Update Beta while A and B are fixed
      VDsq = t(Dsq * t(V))
      UD = t(Beta$d * t(Beta$u))
      xbeta.obs <- suvC(X %*% Beta$v, UD, irow, pcol)
      M_obs = suvC(U, VDsq, irow, pcol)
      S@x = y@x - M_obs - xbeta.obs
      beta = t(X1 %*% S + X2 %*% Beta$v %*% t(UD))
      Beta = fast.svd(as.matrix(beta))
      ##--------------------------------------------
      
      # part 1: Update B
      #VDsq = t(Dsq * t(V))
      IUH = t(U) - (t(U) %*% svdH$u) %*% (svdH$v)
      B = as.matrix(IUH %*% S + (IUH %*% U) %*% t(VDsq))
      if (lambda > 0)
        B = t((B) * (Dsq / (Dsq + lambda)))
      Bsvd = fast.svd(B)
      V = Bsvd$u
      Dsq = Bsvd$d
      U = U %*% (Bsvd$v)
      #--------------------------------
      #-------------------------------------------------------------
      # part 2: Update A
      UDsq = t(Dsq * t(U))
      A.partial = ((S %*% V) + UDsq)
      A = as.matrix(A.partial - svdH$u %*% (svdH$v %*% A.partial))
      if (lambda > 0)
        A = t(t(A) * (Dsq / (Dsq + lambda)))
      Asvd =  fast.svd(A)
      U = Asvd$u
      Dsq = Asvd$d
      V = V %*% (Asvd$v)
      #--------------------------
      # part 3: Update Xbeta
      
      #------------------------------------------------------------
      # part 4: update beta
      # HU = svdH$u %*% (svdH$v %*% U)
      # xbeta.obs = xbeta.obs +
      #   suvC(svdH$u, t(as.matrix(svdH$v %*% S)), irow, pcol) +
      #   suvC(HU, VDsq, irow, pcol)
      
      # VDsq = t(Dsq * t(V))
      # UD = t(Beta$d * t(Beta$u))
      # xbeta.obs <- suvC(X %*% Beta$v, UD, irow, pcol)
      # M_obs = suvC(U, VDsq, irow, pcol)
      # S@x = y@x - M_obs - xbeta.obs
      # beta = t(X1 %*% S + X2 %*% Beta$v %*% t(UD))
      # Beta = fast.svd(as.matrix(beta))
      
      #------------------------------------------------------------------------------
      ratio =  Frob(U.old, Dsq.old, V.old, U, Dsq, V)
      #------------------------------------------------------------------------------
      if (trace.it) {
        obj = (.5 * sum(S@x ^ 2) + lambda * sum(Dsq)) / nz
        if (trace.it)
          cat(iter, ":", "obj", format(round(obj, 5)), "ratio", ratio, "\n")
      }
      #-----------------------------------------------------------------------------------
      
    }
    if (iter == maxit &
        trace.it)
      warning(paste("Convergence not achieved by", maxit, "iterations"))
    if (lambda > 0 & final.svd) {
      #---- update A
      UDsq = t(Dsq * t(U))
      A.partial = ((S %*% V) + UDsq)
      A = as.matrix(A.partial - svdH$u %*% (svdH$v %*% A.partial))
      Asvd =  fast.svd(A)
      U = Asvd$u
      Dsq = Asvd$d
      V = V %*% (Asvd$v)
      #------------------
      Dsq = pmax(Dsq - lambda, 0)
      if (trace.it) {
        UDsq = t(Dsq * t(U))
        M_obs = suvC(UDsq, V, irow, pcol)
        S@x = y@x - M_obs - xbeta.obs
        obj = (.5 * sum(S@x ^ 2) + lambda * sum(Dsq)) / nz
        cat("final SVD:", "obj", format(round(obj, 5)), "\n")
        cat("Number of Iterations for covergence is ", iter, "\n")
      }
    }
    J = min(sum(Dsq > 0) + 1, J)
    J = min(J, length(Dsq))
    out = list(
      u = U[, seq(J), drop = FALSE],
      d = Dsq[seq(J)],
      v = V[, seq(J), drop = FALSE],
      lambda = lambda,
      J = J,
      n_iter = iter,
      xbeta.obs = xbeta.obs,
      beta = beta
    )
    out
  }
