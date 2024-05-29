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
CASMC_fit_rank <-
  function(y,
           X,
           svdH = NULL,
           Xterms = NULL,
           J = 2,
           r = NULL,
           lambda = 0,
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
           init = "naive",
           eps = 1e-17) {
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
      svdH = reduced_hat_decomp(X, 1e-2)
      J_H = svdH$rank
      svdH = svdH$svdH
      if (trace.it)
        print(paste("Rank of H is ", J_H))
    }
    #---------------------------------------------------
    # if Xterms are not provided but X is given.
    if (is.null(Xterms))
      Xterms = GetXterms(X)
    #---------------------------------------------------
    # warm start or initialize (naive or random)
    warm = FALSE
    clean.warm.start(warm.start)
    if (!is.null(warm.start)) {
      #must have u,d and v components
      if (!all(match(c("u", "d", "v", "beta"), names(warm.start), 0) >
               0))
        stop("warm.start does not have components u, d, v, or beta")
      warm = TRUE
      D = warm.start$d
      JD = sum(D > eps)
      beta = warm.start$beta
      if (JD >= J) {
        U = warm.start$u[, seq(J), drop = FALSE]
        V = warm.start$v[, seq(J), drop = FALSE]
        Dsq = D[seq(J)]
      } else{
        Ja = J - JD
        #Dsq = c(D, rep(0, Ja))
        Dsq = c(D, rep(D[JD], Ja))
        U = warm.start$u
        Ua = matrix(rnorm(n * Ja), n, Ja)
        Ua = Ua - U %*% (t(U) %*% Ua)
        Ua = tryCatch(fast.svd(Ua)$u, error = function(e) svd(Ua)$u)
        U = cbind(U, Ua)
        V = cbind(warm.start$v, matrix(0, m, Ja))
        # U_aug = matrix(rnorm(n*Ja), n, Ja)
        # U_aug = qr.Q(qr(cbind(warm.start$u, U_aug)))[,(JD+1):J]
        # V_aug = matrix(rnorm(m*Ja), m, Ja)
        # V_aug = qr.Q(qr(cbind(warm.start$v, V_aug)))[,(JD+1):J]
        # U = cbind(warm.start$u, U_aug)
        # V = cbind(warm.start$v, V_aug)
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
        M <-  tryCatch(propack.svd(M, J),
                       error = function(e){
                         message(paste("Naive Init:",e))
                         svd_trunc_simple(M, J)
                       })
        U = M$u
        V = M$v
        Dsq = pmax(M$d, eps)
        #----------------------
        # initialization for beta = X^-1 Y
        # comment for later: shouldn't be X^-1 H Y??
        beta = as.matrix(ginv(X) %*% Xbeta)
        Y_naive <- Xbeta <- M <- NULL
    print(r)
        #---------------------------------------------------------------
      }
    }
    #----------------------------------------
    yobs <- y@x
    ratio <- 1
    iter <- 0

    if (!is.null(r) && r == 0) {
      beta <- matrix(0, k, m)
      xbeta.obs <- rep(0, length(y@x))
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
      
      if (is.null(r) || r > 0) {
        if (!is.null(r) && r < k && k >= 3) {
          tryCatch({
            beta <- unsvd(svds(beta, r), make_sparse = T)
          }, error = function(e) {
            message(e)
            message(". svds() failed, switching to svd()")
            beta_svd = svd(beta)
            beta_svd$d <- beta_svd$d[1:r, drop=FALSE]
            beta_svd$u <- beta_svd$u[, 1:r, drop=FALSE]
            beta_svd$v <- beta_svd$v[, 1:r, drop=FALSE]
            beta <- unsvd(beta_svd, make_sparse = T)
          })
          
        } else{
          beta = unsvd(tryCatch(fast.svd(beta), error = function(e) svd(beta)))
        }
        if (iter == 1)
          y@x = yobs - M_obs - suvC(X, t(beta), irow, pcol)
        xbeta.obs <- suvC(X, t(beta), irow, pcol)
        beta = as.matrix(beta + Xterms$X1 %*% y)
      }
      
      
      y@x = yobs - M_obs - xbeta.obs
      
      ##--------------------------------------------
      # part 2: Update B while A and beta are fixed
      # updates U, Dsq, V
      B = t(U) %*% y + t(VDsq)
      if (laplace.b)
        B = B - t(V)  %*% L.b
      B = t((B) * (Dsq / (Dsq + lambda)))
      tryCatch({
      Bsvd = fast.svd(as.matrix(B))
      }, error = function(e) message(paste("Bsvd:",e)))
      V = Bsvd$u
      Dsq = pmax(Bsvd$d, eps)
      U = U %*% (Bsvd$v)
      #-------------------------------------------------------------
      # part 3: Update A while B and beta are fixed
      # updates U, Dsq, V
      A = (y %*% V) + t(Dsq * t(U))
      if (laplace.a)
        A = A - L.a %*% U
      A = t(t(A) * (Dsq / (Dsq + lambda)))
      tryCatch({
      Asvd =  fast.svd(as.matrix(A))
      }, error = function(e) message(paste("Asvd:",e)))
      U = Asvd$u
      Dsq = pmax(Asvd$d,eps)
      V = V %*% (Asvd$v)
      #------------------------------------------------------------------------------
      ratio =  Frob(U.old, Dsq.old, V.old, U, Dsq, V)
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
    #print("e")
    # one final fit for one of the parameters (A) has proved to improve the performance significantly.
    if (final.svd) {
      #---- update A
      A = (y %*% V) + t(Dsq * t(U))
      if (laplace.a)
        A = A - L.a %*% U
      A = t(t(A) * (Dsq / (Dsq + lambda)))
      tryCatch({
        Asvd =  fast.svd(as.matrix(A))
      }, error = function(e) message(paste("Asvd:",e)))
      U = Asvd$u
      Dsq = pmax(Asvd$d, eps)
      V = V %*% (Asvd$v)
      Dsq = pmax(Dsq - lambda, eps)
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
    J = min(sum(Dsq > 0) + 1, J)
    J = min(J, length(Dsq))
    if (!is.null(r) && r < k && r > 0 && k >= 3) {
      #beta <- unsvd(svds(beta, r), make_sparse = T)
      tryCatch({
        beta <- unsvd(svds(beta, r), make_sparse = T)
      }, error = function(e) {
        message(e)
        message(". svds() failed, switching to svd()")
        beta_svd = svd(beta)
        beta_svd$d <- beta_svd$d[1:r, drop=FALSE]
        beta_svd$u <- beta_svd$u[, 1:r, drop=FALSE]
        beta_svd$v <- beta_svd$v[, 1:r, drop=FALSE]
        beta <- unsvd(beta_svd, make_sparse = T)
      })
    }
    #print("ee")
    
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
