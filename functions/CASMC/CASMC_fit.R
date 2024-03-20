

CASMC_fit <-
  function (y,
            X = NULL,
            svdH = NULL,
            J = 2,
            thresh = 1e-05,
            lambda = 0,
            maxit = 100,
            trace.it = FALSE,
            warm.start = NULL,
            final.svd = TRUE,
            init = "naive") {
    if (!inherits(y, "dgCMatrix"))
      y = as(y, "dgCMatrix")
    irow = y@i
    pcol = y@p
    n <- dim(y)
    m <- n[2]
    n <- n[1]
    nz = nnzero(y)
    initialize_beta = FALSE
    beta.obs = NULL
    #-------------------------------
    if (is.null(svdH)) {
      stopifnot(!is.null(X))
      svdH = reduced_hat_decomp(X, 1e-2)
      J_H = svdH$rank
      svdH = svdH$svdH
      if (trace.it)
        print(paste("Rank of H is ", J_H))
    }
    #---------------------------------------------------
    warm = FALSE
    clean.warm.start(warm.start)
    if (!is.null(warm.start)) {
      #must have u,d and v components
      if (!all(match(c("u", "d", "v", "xbeta.obs"), names(warm.start), 0) >
               0))
        stop("warm.start does not have components u, d and v")
      warm = TRUE
      D = warm.start$d
      JD = sum(D > 0)
      #J = JD
      xbeta.obs = warm.start$xbeta.obs
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
      stopifnot(init %in% c("random", "naive"))
      if (init == "random") {
        V = matrix(0, m, J)
        U = matrix(rnorm(n * J), n, J)
        U = fast.svd(U)$u
        Dsq = rep(1, J)# we call it Dsq because A=UD and B=VDsq and AB'=U Dsq V^T
        xbeta.obs = suvC(svdH$u, t(as.matrix(svdH$v %*% y)), irow, pcol)
      } else if (init == "naive") {
        Y_naive = as.matrix(y)
        yobs = !is.na(Y_naive)
        Y_naive = naive_MC(Y_naive)
        naive_fit <-  svdH$u %*% (svdH$v  %*% Y_naive)
        xbeta.obs <- naive_fit[yobs]
        naive_fit <- Y_naive - naive_fit
        naive_fit <- propack.svd(as.matrix(naive_fit), J)
        U = naive_fit$u
        V = naive_fit$v
        Dsq = naive_fit$d
      }
    }
    
    #----------------------------------------
    S = y
    ratio <- 1
    iter <- 0
    counter = 0
    best_score = Inf
    best_iter = NA
    #-----------------------------------------------------------
    #####
    # update Beta once at the beginning
    VDsq = t(Dsq * t(V))
    if (is.null(warm.start)) {
      HU = svdH$u %*% (svdH$v %*% U)
      xbeta.obs = xbeta.obs +
        suvC(svdH$u, t(as.matrix(svdH$v %*% S)), irow, pcol) +
        suvC(HU, VDsq, irow, pcol)
    }
    M_obs = suvC(U, VDsq, irow, pcol)
    S@x = y@x - M_obs - xbeta.obs
    #----------------------------------------
    while ((ratio > thresh) & (iter < maxit)) {
      iter <- iter + 1
      U.old = U
      V.old = V
      Dsq.old = Dsq
      #----------------------------------------------
      # part 1: Update B
      VDsq = t(Dsq * t(V))
      UtH = (t(U) %*% svdH$u) %*% (svdH$v)
      IUH = t(U) - UtH
      B = as.matrix(IUH %*% S + (IUH %*% U) %*% t(VDsq))
      if (lambda > 0)
        B = t((B) * (Dsq / (Dsq + lambda)))
      Bsvd = fast.svd(B)
      V = Bsvd$u
      Dsq = Bsvd$d
      U = U %*% (Bsvd$v)
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
      # part 3: Update beta
      VDsq = t(Dsq * t(V))
      HU = svdH$u %*% (svdH$v %*% U)
      
      xbeta.obs = xbeta.obs +
        suvC(svdH$u, t(as.matrix(svdH$v %*% S)), irow, pcol) +
        suvC(HU, VDsq, irow, pcol)
      
      M_obs = suvC(U, VDsq, irow, pcol)
      S@x = y@x - M_obs - xbeta.obs
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
      xbeta.obs = xbeta.obs
    )
    out
  }
