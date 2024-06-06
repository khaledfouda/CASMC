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
CASMC2_fit2 <-
  function(y,
           X,
           svdH = NULL,
           Xterms = NULL,
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
           maxit = 300,
           thresh = 1e-05,
           trace.it = FALSE,
           warm.start = NULL,
           final.svd = TRUE,
           init = "naive",
           min_eigv = 1e-17) {
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
      L.a = computeLaplacian(S.a, normalized = F) * lambda.a
    }
    if (!is.null(S.b) && lambda.b > 0) {
      laplace.b = T
      L.b = computeLaplacian(S.b, normalized = F) * lambda.b
    }
    #--------------------------------
    # if svdH is not given but X is given. Only needed if warm.start is not provided
    if (is.null(warm.start) & is.null(svdH)) {
      svdH = reduced_hat_decomp.H(X)
    }
    #---------------------------------------------------
    # if Xterms are not provided.
    if (is.null(Xterms))
      Xterms = GetXterms(X, lambda.beta)
    #---------------------------------------------------
    # warm start or initialize (naive or random)
    warm = FALSE
    clean.warm.start(warm.start)
    if (!is.null(warm.start)) {
      #must have u,d and v components
      if (!all(match(c("u", "d", "v", "Q", "R"), names(warm.start), 0) >
               0))
        stop("warm.start does not have components u, d, v, or beta")
      warm = TRUE
      # prepare U, Dsq, V
      D = warm.start$d
      JD = sum(D > min_eigv)
      if (JD >= J) {
        U = warm.start$u[, seq(J), drop = FALSE]
        V = warm.start$v[, seq(J), drop = FALSE]
        Dsq = D[seq(J)]
      } else{
        Ja = J - JD
        Dsq = c(D, rep(D[JD], Ja))
        U = warm.start$u
        Ua = matrix(rnorm(n * Ja), n, Ja)
        Ua = Ua - U %*% (t(U) %*% Ua)
        Ua = tryCatch(
          fast.svd(Ua)$u,
          error = function(e)
            svd(Ua)$u
        )
        U = cbind(U, Ua)
        V = cbind(warm.start$v, matrix(0, m, Ja))
      }
      #----------------
      # prepare Q, and R
      Q = warm.start$Q
      R = warm.start$R
      rD = ncol(Q)
      if(rD != r){
        QRsvd = fast.svd(Q%*% t(R))
        if(rD > r){
          QRsvd$d <- diag(sqrt(QRsvd$d[seq(r)]), r, r)
          Q = QRsvd$u[, 1:r, drop=FALSE]  %*% QRsvd$d
          R = QRsvd$v[, 1:r, drop=FALSE] %*% QRsvd$d
        }else{
          ra = r - rD
          Dr = diag(sqrt(c(QRsvd$d, rep(QRsvd$d[ra], ra))), r, r)
          Ur = QRsvd$u
          Ua = matrix(rnorm(k * ra), k, ra)
          Ua = Ua - Ur %*% (t(Ur) %*% Ua)
          Ua = tryCatch(
            fast.svd(Ua)$u,
            error = function(e)
              svd(Ua)$u
          )
          Ur = cbind(Ur, Ua)
          Vr = cbind(QRsvd$v, matrix(0, m, ra))
          Q = Ur  %*% Dr
          R = Vr %*% Dr
        }
      }
      
    } else{
      # initialize. Warm start is not provided
      stopifnot(init %in% c("random", "naive"))
      if (init == "random") {
        stop("Not implemented yet.")
      } else if (init == "naive") {
        Y_naive = as.matrix(y)
        Y_naive = naive_MC(Y_naive)
        Xbeta <-  svdH$u %*% (svdH$v  %*% Y_naive)
        M <- as.matrix(Y_naive - Xbeta)
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
        #----------------------
        # initialization for beta = X^-1 Y
        # comment for later: shouldn't be X^-1 H Y??
        beta = as.matrix(ginv(X) %*% Xbeta)
        QRsvd = fast.svd(beta)
        QRsvd$d <- diag(sqrt(QRsvd$d[seq(r)]), r, r)
        Q = QRsvd$u[, 1:r, drop=FALSE]  %*% sqrt(QRsvd$d)
        R = QRsvd$v[, 1:r, drop=FALSE] %*% sqrt(QRsvd$d)
        Y_naive <- Xbeta <- M <- NULL
        print(r)
        #---------------------------------------------------------------
      }
      Qsvd = fast.svd(Q)
      XtX = t(X) %*% X
    }
    #----------------------------------------
    yobs <- y@x # y will hold the model residuals
    ratio <- 1
    iter <- 0
    print("hi")
    if (!is.null(r) && r == 0) {
      beta <- matrix(0, k, m)
      xbeta.obs <- rep(0, length(y@x))
    }
    print(paste("Q-d", paste(Qsvd$d, collapse = " ")))
    #----------------------------------------
    while ((ratio > thresh) & (iter < maxit)) {
      iter <- iter + 1
      U.old = U
      V.old = V
      Dsq.old = Dsq
      print(paste("b0begin", iter))
      #----------------------------------------------
      # part 0: Update y (training residulas)
      # prereq: U, Dsq, V, Q, R, yobs
      # updates: y; VDsq; XQ
      VDsq = t(Dsq * t(V))
      XQ = X %*% Q
      M_obs = suvC(U, VDsq, irow, pcol)
      xbeta.obs <- suvC(XQ, t(R), irow, pcol)
      y@x = yobs - M_obs - xbeta.obs
      print(paste("b0end", iter))
      #----------------------------------------------
      # part 1: Update R
      # prereq: Q, R, XtX, lambda.beta, y, X, r, Qsvd
      # updates: Q, R
      part1 = t(Q) %*% XtX %*% Q + lambda.beta * diag(1, r, r)
      part2 =  t(XQ) %*% y + (t(Q) %*% (XtX %*% Q)) %*% t(R)
      
      R = as.matrix(ginv(part1) %*% part2)
      R = diag(sqrt(Qsvd$d),r,r) %*% R 
      Rsvd = fast.svd(t(R))
      print(paste("R-d", paste(Rsvd$d, collapse = " ")))
      Q = Qsvd$u %*% Rsvd$v %*% diag(sqrt(Rsvd$d), r, r)
      R = Rsvd$u %*% diag(sqrt(Rsvd$d), r, r)
      print(paste("R1-d", paste(fast.svd(R)$d, collapse = " ")))
      print(paste("Q1-d", paste(fast.svd(Q)$d, collapse = " ")))
      print(paste("p1end", iter))
      #-----------------------------------------------------------------
      # part 2: update Q
      # prereq: Q, R, X, lambda.beta, XtX, y, Rsvd
      # updates: Q, R
      #Q0 = Q
      # for(itt in 1:10)
      Q = (1/ (1+lambda.beta)) * 
        as.matrix(Q + t(X) %*% y %*% R)
      Q = Q %*% diag(sqrt(Rsvd$d),r,r) # QD
      # Q = (1/ (1+lambda.beta)) * 
      #   as.matrix(Q + t(X) %*% y %*% R + XtX %*% (Q0-Q) %*% diag(Rsvd$d^2,r,r))
        #(XtX %*% Q) %*% (t(R) %*% R) 
      print(paste("P2Mid",iter))
      Qsvd = fast.svd(Q)
      print(paste("Q-d", paste(Qsvd$d, collapse = " ")))
      
      Q = Qsvd$u * sqrt(Qsvd$d) #%*% diag(Qsvd$d, r, r)
      R = (Rsvd$v %*% Qsvd$v) * sqrt(Qsvd$d) #diag(Qsvd$d, r, r)
      print(paste("R2-d", paste(fast.svd(R)$d, collapse = " ")))
      print(paste("Q2-d", paste(fast.svd(Q)$d, collapse = " ")))
      #-------------------------------------------------------------------
      # part extra: re-update y 
      # xbeta.obs <- suvC(as.matrix(X %*% Q), as.matrix(t(R)), irow, pcol)
      # y@x = yobs - M_obs - xbeta.obs
      print(paste("p2end", iter))
      #--------------------------------------------------------------------
      ##--------------------------------------------
      # part 3: Update B 
      # prereq: U, VDsq, y, Dsq, lambda.M, L.b
      # updates U, Dsq, V
      B = t(U) %*% y + t(VDsq)
      if (laplace.b)
        B = B - t(V)  %*% L.b
      B = as.matrix(t((B) * (Dsq / (Dsq + lambda.M))))
      Bsvd = tryCatch({
        fast.svd(B)
      }, error = function(e){
        message(paste("Loop/B:", e))
        svd(B)
      })
      V = Bsvd$u
      Dsq = pmax(Bsvd$d, min_eigv)
      U = U %*% (Bsvd$v)
      print(paste("p3end", iter))
      #-------------------------------------------------------------
      # part 4: Update A 
      # prereq: U, D, VDsq, y, lambda.M, L.a
      # updates U, Dsq, V
      A = (y %*% V) + t(Dsq * t(U))
      if (laplace.a)
        A = A - L.a %*% U
      A = as.matrix(t(t(A) * (Dsq / (Dsq + lambda.M))))
      Asvd = tryCatch({
        fast.svd(A)
      }, error = function(e){
        message(paste("Loop/A:", e))
        svd(A)
      })
      U = Asvd$u
      Dsq = pmax(Asvd$d, min_eigv)
      V = V %*% (Asvd$v)
      print(paste("p4end", iter))
      #------------------------------------------------------------------------------
      ratio =  Frob(U.old, Dsq.old, V.old, U, Dsq, V)
      print(paste("ratio:", ratio))
      #------------------------------------------------------------------------------
      if (trace.it) {
        obj = (.5 * sum(y@x ^ 2) + lambda.M * sum(Dsq)) / nz
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
      A = as.matrix(t(t(A) * (Dsq / (Dsq + lambda.M))))
      Asvd = tryCatch({
        fast.svd(A)
      }, error = function(e){
        message(paste("Final/A:", e))
        svd(A)
      })
      U = Asvd$u
      V = V %*% Asvd$v
      Dsq = pmax(Asvd$d - lambda.M, min_eigv)
      #------------------
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
    J = min(sum(Dsq > 0) + 1, J)
    J = min(J, length(Dsq))
    
    out = list(
      u = U[, seq(J), drop = FALSE],
      d = Dsq[seq(J)],
      v = V[, seq(J), drop = FALSE],
      Q = Q,
      R = R,
      lambda.M = lambda.M,
      J = J,
      lambda.beta = lambda.beta,
      lambda.a = lambda.a,
      lambda.b = lambda.b,
      n_iter = iter
    )
    out
  }

