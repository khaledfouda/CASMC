






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
  CASMC_Nuclear_fit <-
  function(y,
           X,
           XtX = t(X) %*% X,
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
    laplace.a = laplace.b = F
    if (!is.null(S.a) && lambda.a > 0) {
      laplace.a = T
      L.a = utils$computeLaplacian(S.a, normalized = normalized_laplacian) * lambda.a
    }
    if (!is.null(S.b) && lambda.b > 0) {
      laplace.b = T
      L.b = utils$computeLaplacian(S.b, normalized = normalized_laplacian) * lambda.b
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
      warm.out <-
        utils$prepare.M.warm.start(warm.start, J, n, m, min_eigv)
      U = warm.out$U
      V = warm.out$V
      Dsq = warm.out$Dsq
      #----------------
      # prepare ub, db, vb
      warm.out <-
        utils$prepare.M.warm.start(list(
          u = warm.start$ub,
          v = warm.start$vb,
          d = warm.start$db
        ),
        r,
        k,
        m,
        min_eigv)
      Ub = warm.out$U
      Vb = warm.out$V
      Db = warm.out$Dsq
      # print(paste(r, dim(Ub)[2], dim(warm.start$ub)[2],length(Db),
      #             length(warm.start$db)))
      Q = utils$UD(Ub, Db)
      R = utils$UD(Vb, Db)
    } else{
      #---------------------------------------------------------------
      Y_naive = naive_MC(as.matrix(y))
      beta = utils$inv(XtX, T) %*% t(X) %*% Y_naive
      QRsvd = utils$svdopt(beta, r, k, m, T, F)
      Ub <- as.matrix(QRsvd$u, k, r)
      Db <- pmax(sqrt(QRsvd$d), min_eigv)
      Vb <- as.matrix(QRsvd$v, m, r)
      Q = utils$UD(Ub, Db)
      R = utils$UD(Vb, Db)
      M <- as.matrix(Y_naive - (X %*% Q) %*% t(R))
      M <- utils$svdopt(M, J, n, m, F, F)
      U = M$u
      V = M$v
      Dsq = M$d
      Y_naive <- beta <- M <- NULL
      #----------------------
      #---------------------------------------------------------------
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
      RD = t(as.matrix(utils$inv(part1, T) %*% part2) * Db)
      Rsvd = utils$svd_small_nc(RD, F, p = r)
      Ub = Ub %*% Rsvd$v
      Db <- sqrt(Rsvd$d) #pmax(sqrt(Rsvd$d), min_eigv)
      Vb = Rsvd$u
      Q = utils$UD(Ub , Db)
      R = utils$UD(Vb , Db)
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
      part1 <-
        as.vector(t(X) %*% y %*% R + XtX %*% utils$UD(Q, Db ^ 2))
      part2 <-
        kronecker(diag(Db ^ 2, r, r), XtX) + diag(lambda.beta, k * r, k * r)
      Q <- matrix(utils$inv(part2, T) %*% part1, k, r)
      
      Qsvd = utils$svd_small_nc(as.matrix(utils$UD(Q, Db)), F, p = r)
      # Qsvd = utils$fast.svd(as.matrix(utils$UD(Q, Db)), trim = FALSE)
      Ub = Qsvd$u
      Db <- sqrt(Qsvd$d) #pmax(sqrt(Qsvd$d), min_eigv)
      Vb = Vb %*% Qsvd$v
      Q = utils$UD(Ub, Db)
      R = utils$UD(Vb, Db)
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
      Bsvd = utils$svd_small_nc(B, FALSE, p = J) 
      V = Bsvd$u
      Dsq = Bsvd$d #pmax(Bsvd$d, min_eigv)
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
      Asvd =  utils$svd_small_nc(A, FALSE, p = J)
      U = Asvd$u
      Dsq = Asvd$d #pmax(Asvd$d, min_eigv)
      V = V %*% (Asvd$v)
      #-------------------------------------------------------------------
      # part extra: re-update y
      # VDsq = t(Dsq * t(V))
      # M_obs = suvC(U, VDsq, irow, pcol)
      # y@x = ypart - M_obs
      # #------------------------------------------------------------------------------
      ratio =  utils$Frob(U.old, Dsq.old, V.old, U, Dsq, V) +
        utils$Frob(Ub.old, Db.old, Vb.old, Ub, Db, Vb)
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
      Asvd =  utils$svd_small_nc(A, FALSE, p = J)
      U = Asvd$u
      V = V %*% Asvd$v
      Dsq = pmax(Asvd$d - lambda.M, 0)
      
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
    J = min(max(1, sum(Dsq > min_eigv)), J)
    r = min(max(1, sum(Db > min_eigv)), r)
    
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

