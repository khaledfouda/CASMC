

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
           Qtype = 1,
           qiter.max = 10,
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
    # warm start or initialize (naive or random)
    warm = FALSE
    clean.warm.start(warm.start)
    if (!is.null(warm.start)) {
      #must have u,d and v components
      if (!all(match(c("u", "d", "v", "ub", "db", "vb"), names(warm.start), 0) >
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
      Db = warm.start$db
      rD = sum(diag(Db) > 0)
      if (rD >= r) {
        Ub = warm.start$ub[, seq(r), drop = FALSE]
        Vb = warm.start$vb[, seq(r), drop = FALSE]
        Db = Db[seq(r),seq(r), drop = FALSE]
      } else{
        ra = r - rD
        Db = diag(c(diag(Db), rep(diag(Db)[rD], ra)), r, r)
        Ub = warm.start$ub
        Uba = matrix(rnorm(k * ra), k, ra)
        Uba = Uba - Ub %*% (t(Ub) %*% Uba)
        Uba = tryCatch(
          fast.svd(Uba)$u,
          error = function(e)
            svd(Uba)$u
        )
        Ub = cbind(Ub, Uba)
        Vb = cbind(warm.start$vb, matrix(0, m, ra))
        # print(dim(Ub))
        # print(dim(Vb))
        # print(dim(Db))
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
        QRsvd = svd_trunc_simple(beta, r)
        Ub <- QRsvd$u
        Db <- diag(sqrt(QRsvd$d), r, r)
        Vb <- QRsvd$v
        Y_naive <- Xbeta <- M <- NULL
        #---------------------------------------------------------------
      }
    }
    Q = Ub %*% Db
    R = Vb %*% Db
    XtX = t(X) %*% X
    #----------------------------------------
    yobs <- y@x # y will hold the model residuals
    ratio <- 1
    iter <- 0
    if (!is.null(r) && r == 0) {
      beta <- matrix(0, k, m)
      xbeta.obs <- rep(0, length(y@x))
    }
    #print(paste("Q-d", paste(Qsvd$d, collapse = " ")))
    #----------------------------------------
    while ((ratio > thresh) & (iter < maxit)) {
      iter <- iter + 1
      U.old = U
      V.old = V
      Dsq.old = Dsq
      # print(paste("b0begin", iter))
      #----------------------------------------------
      # part 0: Update y (training residulas)
      # prereq: U, Dsq, V, Q, R, yobs
      # updates: y; VDsq; XQ
      VDsq = t(Dsq * t(V))
      XQ = X %*% Q
      M_obs = suvC(U, VDsq, irow, pcol)
      if(iter == 1){
      xbeta.obs <- suvC(XQ, R, irow, pcol)
      y@x = yobs - M_obs - xbeta.obs}
      # print(paste("b0end", iter))
      #----------------------------------------------
      # part 1: Update R
      # prereq: Q, R, XtX, lambda.beta, y, X, r
      # updates: Q, R
      part1 = t(Q) %*% XtX %*% Q + lambda.beta * diag(1, r, r)
      part2 =  t(XQ) %*% y + (t(Q) %*% (XtX %*% Q)) %*% t(R)
      
      RD = t(as.matrix(ginv(part1) %*% part2)) %*% Db
      Rsvd = svd(RD)
      Ub = Ub %*% Rsvd$v
      # print(dim(Db))
      # print(length(Rsvd$d))
      Db[cbind(1:r,1:r)] <- sqrt(Rsvd$d)
      Vb = Rsvd$u
      
      Q = Ub %*% Db 
      R = Vb %*% Db
      
      
      # print(paste("R-d", paste(Rsvd$d, collapse = " ")))
      # # Q = Qsvd$u %*% Rsvd$v %*% diag(sqrt(Rsvd$d), r, r)
      # # R = Rsvd$u %*% diag(sqrt(Rsvd$d), r, r)
      # print(paste("R1-d", paste(fast.svd(R)$d, collapse = " ")))
      # print(paste("Q1-d", paste(fast.svd(Q)$d, collapse = " ")))
      # print(paste("p1end", iter))
      #-----------------------------------------------------------------
      # part 2: update Q
      # prereq: Q, R, X, lambda.beta, XtX, y, Rsvd
      # updates: Q, R
      # Q0 = Qold = Q
      # for(itt in 1:60){
      # #y@x = yobs
      # #gradient = t(X) %*% y %*% R - XtX %*% Q %*% t(R) %*% R - t(X) %*% U %*% t(VDsq) %*% R
      # 
      # #Q = (1/ (1+0.1*lambda.beta)) * 
      # #  as.matrix(Q + 0.1*gradient)
      #   # as.matrix(Q + t(X) %*% y %*% R)
      # 
      # #Q = Q %*% diag(sqrt(Rsvd$d),r,r) # QD
      # #print(paste(paste(dim(R),collapse=" "),paste(dim(Q),collapse=" ")))
      # 
      #    Q = (1/ (1+lambda.beta)) * 
      #    as.matrix(Q + t(X) %*% y %*% R + XtX %*% (Q0-Q) %*% (Db^2) )
      #    print(round(sum( (Q-Qold)^2 ),3))
      #    Qold  = Q
      # }
      Q0 =  Q
      Qold = matrix(0, nrow(Q), ncol(Q))
      qiter = 0
      Dbsq = Db^2
      # newton!!
      hessian <- matrix(0,k,r)
      while(sqrt(sum((Q - Qold)^2)) > 1e-3 & qiter < qiter.max ){
        Qold <- Q
        if(Qtype == 1){
          
        gradient = - t(X) %*% y %*% R + XtX %*% (Q-Q0) %*% (Dbsq) + lambda.beta * Q
        for(ih in 1:k)
          hessian[ih,] <- sum(XtX[ih,]) * diag(Dbsq)
        hessian = hessian + diag(lambda.beta, k, r)

        Q <- Q -  gradient / hessian
        }else{
          # proximal descent
          Q = (1/ (1+lambda.beta)) * 
                as.matrix(Q + t(X) %*% y %*% R + XtX %*% (Q0-Q) %*% (Dbsq) )
        }
        qiter <- qiter + 1
      }
      Qsvd = svd(Q %*% Db)
      
      # print(dim(Db))
      # print(length(Qsvd$d))
      
      Ub = Qsvd$u
      Db[cbind(1:r,1:r)] <- sqrt(Qsvd$d)
      Vb = Vb %*% Qsvd$v
      
      Q = Ub %*% Db 
      R = Vb %*% Db
      #-------------------------------------------------------------------
      # part extra: re-update y 
      xbeta.obs <- suvC(as.matrix(X %*% Q), as.matrix((R)), irow, pcol)
      y@x = yobs - M_obs - xbeta.obs
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
        svd(B)
      }, error = function(e){
        message(paste("Loop/B:", e))
        svd(B)
      })
      V = Bsvd$u
      Dsq = pmax(Bsvd$d, min_eigv)
      U = U %*% (Bsvd$v)
      #print(paste("p3end", iter))
      #-------------------------------------------------------------
      # part 4: Update A 
      # prereq: U, D, VDsq, y, lambda.M, L.a
      # updates U, Dsq, V
      A = (y %*% V) + t(Dsq * t(U))
      if (laplace.a)
        A = A - L.a %*% U
      A = as.matrix(t(t(A) * (Dsq / (Dsq + lambda.M))))
      Asvd = tryCatch({
        svd(A)
      }, error = function(e){
        message(paste("Loop/A:", e))
        svd(A)
      })
      U = Asvd$u
      Dsq = pmax(Asvd$d, min_eigv)
      V = V %*% (Asvd$v)
      #print(paste("p4end", iter))
      #------------------------------------------------------------------------------
      ratio =  Frob(U.old, Dsq.old, V.old, U, Dsq, V)
      #------------------------------------------------------------------------------
      if (trace.it) {
        obj = (.5 * sum(y@x ^ 2) + lambda.M * sum(Dsq)) / nz
        cat(iter, ":", "obj", format(round(obj, 5)), "ratio", ratio, "qiter",qiter, "\n")
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
      #---------------------------------------------
      # part1 = t(Q) %*% XtX %*% Q + lambda.beta * diag(1, r, r)
      # part2 =  t(XQ) %*% y + (t(Q) %*% (XtX %*% Q)) %*% t(R)
      # 
      # RD = t(as.matrix(ginv(part1) %*% part2)) %*% Db
      # Rsvd = fast.svd(RD)
      # Ub = Ub %*% Rsvd$v
      # Db[cbind(1:r,1:r)] <- sqrt(Rsvd$d)
      # Vb = Rsvd$u
      
      #Db[cbind(1:r,1:r)] <- pmax(Db[cbind(1:r,1:r)] - lambda.beta, 0)
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
    r = min(sum(Db>0) + 1, r)
    r = min(r, length(Db))
    
    out = list(
      u = U[, seq(J), drop = FALSE],
      d = Dsq[seq(J)],
      v = V[, seq(J), drop = FALSE],
      ub = Ub[, seq(r), drop = FALSE],
      db = Db[seq(r),seq(r), drop=FALSE],
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

