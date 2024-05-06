



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
CASMC_fit_laplace <-
  function(y,
           X,
           svdH = NULL,
           Xterms = NULL,
           J = 2,
           r = NULL,
           lambda = 0,
           lambda.a = 0,
           lambda.b = 0,
           S.a = NULL,
           # similarity matrix for A
           S.b = NULL,
           # similarity matrix for B
           maxit = 100,
           thresh = 1e-05,
           trace.it = FALSE,
           warm.start = NULL,
           final.svd = TRUE,
           init = "naive") {
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
    laplace.a= laplace.b = F
    if (!is.null(S.a) & lambda.a > 0) {
      laplace.a = T
      L.a = computeLaplacian(S.a, normalized = TRUE) * lambda.a
    }
    if (!is.null(S.b) & lambda.b > 0) {
      laplace.b = T
      L.b = computeLaplacian(S.b, normalized = TRUE) * lambda.b
    }
    #--------------------------------
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
    Xterms <- NULL
    #---------------------------------------------------
    # warm start or initialize (naive or random)
    warm = FALSE
    clean.warm.start(warm.start)
    if (!is.null(warm.start)) {
      #must have u,d and v components
      if (!all(match(c("u", "d", "v", "xbeta.obs", "Beta"), names(warm.start), 0) >
               0))
        stop("warm.start does not have components u, d and v")
      warm = TRUE
      D = warm.start$d
      JD = sum(D > 0)
      #J = JD
      #beta = as.matrix(warm.start$beta)
      Beta = warm.start$Beta #fast.svd(beta)
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
        # initialization for beta = X^-1 Y
        # comment for later: shouldn't be X^-1 H Y??
        beta = t(ginv(X) %*% Y_naive)
        Beta = fast.svd(beta)
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
    # Adjst the rank of Beta if provided
    if (!is.null(r)) {
      beta_rank = min(sum(round(Beta$d, 4) > 0), r)
      if (beta_rank == 0) {
        Beta$u <- matrix(0, m, 1)
        Beta$v <- matrix(0, k, 1)
        Beta$d <- c(0)
      } else{
        Beta$u <- Beta$u[, 1:beta_rank, drop = FALSE]
        Beta$v <- Beta$v[, 1:beta_rank, drop = FALSE]
        Beta$d <- Beta$d[1:beta_rank]
      }
    }
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
      if (iter == 1) {
        UD.beta = t(Beta$d * t(Beta$u))
        # print(dim(X)); print(dim(Beta$v))
        xbeta.obs <- suvC(X %*% Beta$v, UD.beta, irow, pcol)
        S@x = y@x - M_obs - xbeta.obs
      }
      # print(dim(X1 %*% S))
      # print(dim(X2 %*% Beta$v %*% t(UD.beta)))
      # print(dim(X2))
      # print(dim(t(UD.beta)))
      # print(dim(Beta$v))
      beta = t(X1 %*% S + X2 %*% Beta$v %*% t(UD.beta))
      Beta = fast.svd(as.matrix(beta))
      # Adjust the rank of Beta if provided
      if (!is.null(r)) {
        beta_rank = min(sum(round(Beta$d, 4) > 0), r)
        if (beta_rank == 0) {
          Beta$u <- matrix(0, m, 1)
          Beta$v <- matrix(0, k, 1)
          Beta$d <- c(0)
        } else{
          Beta$u <- Beta$u[, 1:beta_rank, drop = FALSE]
          Beta$v <- Beta$v[, 1:beta_rank, drop = FALSE]
          Beta$d <- Beta$d[1:beta_rank]
          
        }
      }
      # why not this? [added on Mai 2nd - not tested yet; also lines 155 to 163]
      UD.beta = t(Beta$d * t(Beta$u))
      xbeta.obs <- suvC(X %*% Beta$v, UD.beta, irow, pcol)
      S@x = y@x - M_obs - xbeta.obs
      ##--------------------------------------------
      # part 2: Update B while A and beta are fixed
      # updates U, Dsq, V
      #IUH = t(U) - (t(U) %*% svdH$u) %*% (svdH$v)
      B = as.matrix(t(U) %*% S + t(VDsq))
      if(laplace.b) B = B - t(V)  %*% L.b
      #B = as.matrix(IUH %*% S + (IUH %*% U) %*% t(VDsq))
      B = t((B) * (Dsq / (Dsq + lambda)))
      Bsvd = fast.svd(B)
      V = Bsvd$u
      Dsq = Bsvd$d
      U = U %*% (Bsvd$v)
      #-------------------------------------------------------------
      # part 3: Update A while B and beta are fixed
      # updates U, Dsq, V
      #UDsq = t(Dsq * t(U))
      #print(hi)
      A = as.matrix((S %*% V) + t(Dsq * t(U)))
      #A = as.matrix(A.partial - svdH$u %*% (svdH$v %*% A.partial))
      if(laplace.a) A = A - L.a %*% U 
      A = t(t(A) * (Dsq / (Dsq + lambda)))
      Asvd =  fast.svd(A)
      U = Asvd$u
      Dsq = Asvd$d
      V = V %*% (Asvd$v)
      #------------------------------------------------------------------------------
      ratio =  Frob(U.old, Dsq.old, V.old, U, Dsq, V)
      #------------------------------------------------------------------------------
      if (trace.it) {
        obj = (.5 * sum(S@x ^ 2) + lambda * sum(Dsq)) / nz
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
      A.partial = ((S %*% V) + t(Dsq * t(U)))
      A = as.matrix(A.partial - svdH$u %*% (svdH$v %*% A.partial))
      Asvd =  fast.svd(A)
      U = Asvd$u
      Dsq = Asvd$d
      V = V %*% (Asvd$v)
      Dsq = pmax(Dsq - lambda, 0)
      #------------------
      if (trace.it) {
        M_obs = suvC(t(Dsq * t(U)), V, irow, pcol)
        S@x = y@x - M_obs - xbeta.obs
        obj = (.5 * sum(S@x ^ 2) + lambda * sum(Dsq)) / nz
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
      lambda = lambda,
      J = J,
      n_iter = iter,
      xbeta.obs = xbeta.obs,
      Beta = Beta
    )
    out
  }
