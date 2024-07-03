





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
CASMC2_fit <-
  function(y,
           X,
           J = 2,
           r = 2,
           lambda.M = 0,
           lambda.beta = 0,
           # similarity matrix for A
           S.a = NULL,
           lambda.a = 0,
           # similarity matrix for B
           S.b = NULL,
           lambda.b = 0,
           normalized_laplacian = TRUE,
           # convergence parameters
           maxit = 300,
           thresh = 1e-05,
           trace.it = FALSE,
           warm.start = NULL,
           # the following should not be modified
           final.svd = TRUE,
           init = "naive",
           min_eigv = 1e-4) {
    stopifnot(inherits(y, "dgCMatrix"))
    irow = y@i
    pcol = y@p
    n <- dim(y)
    m <- n[2]
    n <- n[1]
    k <- ncol(X)
    XtX = t(X) %*% X
    
    if (trace.it)
      nz = nnzero(y, na.counted = TRUE)
    #-------------------------------
    laplace.a = laplace.b = F
    if (!is.null(S.a) && lambda.a > 0) {
      laplace.a = T
      L.a = computeLaplacian(S.a, normalized = normalized_laplacian) * lambda.a
    }
    if (!is.null(S.b) && lambda.b > 0) {
      laplace.b = T
      L.b = computeLaplacian(S.b, normalized = normalized_laplacian) * lambda.b
    }
    #---------------------------------------------------
    # warm start or initialize (naive or random)
    warm = FALSE
    clean.warm.start(warm.start)
    if (!is.null(warm.start)) {
      #must have u,d and v components
      if (!all(match(c("u", "d", "v", "ub", "db", "vb"), names(warm.start), 0) >
               0))
        stop("warm.start is missing some or all of the components.")
      warm = TRUE
      # prepare U, Dsq, V
      D = warm.start$d
      JD = sum(D > min_eigv)
      if (JD >= J) {
        U = warm.start$u[, seq(J), drop = FALSE]
        V = warm.start$v[, seq(J), drop = FALSE]
        Dsq = D[seq(J)]
      } else{
        # upscale
        Ja = J - JD
        Dsq = c(D, rep(D[JD], Ja))
        U = warm.start$u
        Ua = matrix(rnorm(n * Ja), n, Ja)
        Ua = Ua - U %*% (t(U) %*% Ua)
        Ua = tryCatch(
          fast.svd(Ua, trim = FALSE)$u,
          error = function(e)
            svd(Ua)$u
        )
        U = cbind(U, Ua)
        V = cbind(warm.start$v, matrix(0, m, Ja))
      }
      #----------------
      # prepare ub, db, vb
      Db = warm.start$db
      rD = sum(Db > min_eigv)
      if (rD >= r) {
        Ub = warm.start$ub[, seq(r), drop = FALSE]
        Vb = warm.start$vb[, seq(r), drop = FALSE]
        Db = Db[seq(r)]
      } else{
        ra = r - rD
        Db = c(Db, rep(Db[rD], ra))
        Ub = warm.start$ub
        Uba = matrix(rnorm(k * ra), k, ra)
        Uba = Uba - Ub %*% (t(Ub) %*% Uba)
        Uba = tryCatch(
          fast.svd(Uba, trim = FALSE)$u,
          error = function(e)
            svd(Uba)$u
        )
        Ub = cbind(Ub, Uba)
        Vb = cbind(warm.start$vb, matrix(0, m, ra))
      }
      
      Q = UD(Ub, Db)
      R = UD(Vb, Db)
    } else{
      # initialize. Warm start is not provided
      stopifnot(init %in% c("random", "naive"))
      if (init == "random") {
        stop("Not implemented yet.")
      } else if (init == "naive") {
        # Y_naive = as.matrix(y)
        # Y_naive = naive_MC(Y_naive)
        # svdH = reduced_hat_decomp.H(X)
        # Xbeta <-  svdH$u %*% (svdH$v  %*% Y_naive)
        # M <- as.matrix(Y_naive - Xbeta)
        # M <-  tryCatch(
        #   propack.svd(M, J),
        #   error = function(e) {
        #     message(paste("Naive Init:", e))
        #     svd_trunc_simple(M, J)
        #   }
        # )
        # U = M$u
        # V = M$v
        # Dsq = pmax(M$d, min_eigv)
        # #----------------------
        # # initialization for beta = X^-1 Y
        # # comment for later: shouldn't be X^-1 H Y??
        # beta = as.matrix(ginv(X) %*% Xbeta)
        # QRsvd = svd_trunc_simple(beta, r)
        # Ub <- as.matrix(QRsvd$u, k, r)
        # Db <- pmax(sqrt(QRsvd$d), min_eigv)
        # Vb <- as.matrix(QRsvd$v, m, r)
        # Y_naive <- Xbeta <- M <- NULL
        #---------------------------------------------------------------
        Y_naive = naive_MC(as.matrix(y))
        beta = ginv(XtX) %*% t(X) %*% Y_naive
        QRsvd = svd_trunc_simple(beta, r)
        Ub <- as.matrix(QRsvd$u, k, r)
        Db <- pmax(sqrt(QRsvd$d), min_eigv)
        Vb <- as.matrix(QRsvd$v, m, r)
        Q = UD(Ub, Db)
        R = UD(Vb, Db)
        M <- as.matrix(Y_naive - (X %*% Q) %*% t(R) )
        M <-  tryCatch(
          propack.svd(M, J),
          error = function(e) {
            message(paste("Naive Init:", e))
            svd_trunc_simple(M, J)
          }
        )
        U = M$u
        V = M$v
        Dsq = pmax(M$d, min_eigv)
        Y_naive <- beta <- M <- NULL
        #----------------------
        #---------------------------------------------------------------
      }
    }
    
    #----------------------------------------
    yobs <- y@x # y will hold the model residuals
    ratio <- 1
    iter <- 0
    XQ = X %*% Q
    xbeta.obs <- suvC(XQ, R, irow, pcol)
    ypart = yobs - xbeta.obs
    # if (!is.null(r) && r == 0) { # not used but should be
    #   beta <- matrix(0, k, m)
    #   xbeta.obs <- rep(0, length(y@x))
    # }
    #----------------------------------------
    while ((ratio > thresh) & (iter < maxit)) {
      iter <- iter + 1
      U.old = U
      V.old = V
      Dsq.old = Dsq
      Ub.old = Ub
      Vb.old = Vb
      Db.old = Db
      #----------------------------------------------
      # part 0: Update y (training residulas)
      # prereq: U, Dsq, V, Q, R, yobs
      # updates: y; VDsq; XQ
      VDsq = t(Dsq * t(V))
      # ypart should be  yobs - xbeta.obs at the beginning:
      M_obs = suvC(U, VDsq, irow, pcol)
      y@x = ypart - M_obs
      ypart = yobs - M_obs
      #----------------------------------------------
      # part 1: Update R
      # prereq: Q, R, XtX, lambda.beta, y, X, r
      # updates: Q, R
      QXtXQ <- t(Q) %*% (XtX %*% Q)
      
      part1 = (QXtXQ + diag(lambda.beta, r, r))
      part2 =  (t(XQ) %*% y + QXtXQ %*% t(R))
      RD = t( as.matrix(solve(part1) %*% part2) * Db)
      
      Rsvd = fast.svd(RD, trim = FALSE)
      Ub = Ub %*% Rsvd$v
      Db <- pmax(sqrt(Rsvd$d), min_eigv)
      Vb = Rsvd$u
      Q = UD(Ub , Db)
      R = UD(Vb , Db)
      #-------------------------------------------------------------------
      # part extra: re-update y
      xbeta.obs <-
        suvC(X %*% Q, R, irow, pcol)
      y@x = ypart - xbeta.obs
      #--------------------------------------------------------------------
      #-----------------------------------------------------------------
      # part 2: update Q
      # prereq: Q, R, X, lambda.beta, XtX, y, Rsvd
      # updates: Q, R
      part1 <- as.vector(t(X) %*% y %*% R + XtX %*% UD(Q, Db^2))
      part2 <- kronecker(diag(Db^2, r, r), XtX) + diag(lambda.beta, k*r,k*r)
      Q <- matrix(solve(part2) %*% part1, k, r)
      
      Qsvd = fast.svd(as.matrix(UD(Q, Db)), trim = FALSE)
      Ub = Qsvd$u
      Db <- pmax(sqrt(Qsvd$d), min_eigv)
      Vb = Vb %*% Qsvd$v
      Q = UD(Ub,Db)
      R = UD(Vb,Db)
      XQ = X %*% Q
      #-------------------------------------------------------------------
      # part extra: re-update y
      xbeta.obs <-
        suvC(X %*% Q, R, irow, pcol)
      y@x = ypart - xbeta.obs
      ypart = yobs - xbeta.obs
      #--------------------------------------------------------------------
      # print("hi2")
      ##--------------------------------------------
      # part 3: Update B
      # prereq: U, VDsq, y, Dsq, lambda.M, L.b
      # updates U, Dsq, V
      B = t(U) %*% y + t(VDsq)
      if (laplace.b)
        B = B - t(V)  %*% L.b
      B = as.matrix(t((B) * (Dsq / (Dsq + lambda.M))))
      Bsvd = tryCatch({
        fast.svd(B, trim = FALSE)
      }, error = function(e) {
        message(paste("Loop/B:", e))
        svd(B)
      })
      V = Bsvd$u
      Dsq = pmax(Bsvd$d, min_eigv)
      U = U %*% (Bsvd$v)
      #-------------------------------------------------------------------
      # part extra: re-update y
      VDsq = t(Dsq * t(V))
      M_obs = suvC(U, VDsq, irow, pcol)
      y@x = ypart - M_obs
      #--------------------------------------------------------------------
      #-------------------------------------------------------------
      # part 4: Update A
      # prereq: U, D, VDsq, y, lambda.M, L.a
      # updates U, Dsq, V
      A = (y %*% V) + t(Dsq * t(U))
      if (laplace.a)
        A = A - L.a %*% U
      A = as.matrix(t(t(A) * (Dsq / (Dsq + lambda.M))))
      Asvd = tryCatch({
        #svd(A)
        fast.svd(A, trim = FALSE)
      }, error = function(e) {
        message(paste("Loop/A:", e))
        svd(A)
      })
      U = Asvd$u
      Dsq = pmax(Asvd$d, min_eigv)
      V = V %*% (Asvd$v)
      #-------------------------------------------------------------------
      # part extra: re-update y
      # VDsq = t(Dsq * t(V))
      # M_obs = suvC(U, VDsq, irow, pcol)
      # y@x = ypart - M_obs
      # #------------------------------------------------------------------------------
      ratio =  Frob(U.old, Dsq.old, V.old, U, Dsq, V) +
               Frob(Ub.old, Db.old, Vb.old, Ub, Db, Vb)
      #------------------------------------------------------------------------------
      if (trace.it) {
        obj = (.5 * sum(y@x ^ 2) +
                 lambda.M * sum(Dsq) +
                 lambda.beta * sum(Db ^ 2)) / nz
        cat(iter,
            ":",
            "obj",
            format(round(obj, 5)),
            "ratio",
            ratio,
            "\n")
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
      A = as.matrix(t(t(A) * (Dsq / (Dsq + lambda.M))))
      Asvd = tryCatch({
        fast.svd(A, trim = FALSE)
      }, error = function(e) {
        message(paste("Final/A:", e))
        svd(A)
      })
      U = Asvd$u
      V = V %*% Asvd$v
      # is this good?
      Dsq = pmax(Asvd$d - lambda.M, min_eigv)
      #---------------------------------------------
      # do it again with Q
      # M_obs = suvC(t(Dsq * t(U)), V, irow, pcol)
      # y@x = yobs - M_obs - xbeta.obs
      # 
      # part1 <- as.vector(t(X) %*% y %*% R + XtX %*% UD(Q, Db^2))
      # part2 <- kronecker(diag(Db^2, r, r), XtX) + diag(lambda.beta, k*r,k*r)
      # Q <- matrix(solve(part2) %*% part1, k, r)
      # 
      # Qsvd = fast.svd(as.matrix(UD(Q, Db)), trim = FALSE)
      # Ub = Qsvd$u
      # Db <- sqrt(pmax(Qsvd$d - lambda.beta, min_eigv))
      # Vb = Vb %*% Qsvd$v
      #---------------------------------------------
      if (trace.it) {
        M_obs = suvC(t(Dsq * t(U)), V, irow, pcol)
        y@x = yobs - M_obs - xbeta.obs
        obj = (.5 * sum(y@x ^ 2) + lambda.M * sum(Dsq)) / nz
        cat("final SVD:", "obj", format(round(obj, 5)), "\n")
        cat("Number of Iterations for covergence is ", iter, "\n")
      }
    }
    
    #-------------------------------------------------------
    # trim in case we reduce the rank of M to be smaller than J.
    J = min(max(1,sum(Dsq > min_eigv)), J)
    r = min(max(1,sum(Db > min_eigv)), r)
    
    out = list(
      u = U[, seq(J), drop = FALSE],
      d = Dsq[seq(J)],
      v = V[, seq(J), drop = FALSE],
      ub = Ub[, seq(r), drop = FALSE],
      db = Db[seq(r)],
      vb = Vb[, seq(r), drop = FALSE],
      #---------------------
      # hyperparameters
      lambda.M = lambda.M,
      lambda.beta = lambda.beta,
      J = J,
      r = r,
      lambda.a = lambda.a,
      lambda.b = lambda.b,
      #--------------------------
      # convergence
      n_iter = iter
    )
    out
  }

