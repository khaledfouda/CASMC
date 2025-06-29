




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
#'  CAMC_fit(y,X,J=5)
#' @export
#'
CAMC3_fit <-
  CAMC_Lasso_fit <-
  function(y,
           X,
           XtX = t(X) %*% X,
           J = 2,
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
           beta.iter.max = 10,
           learning.rate = 0.001,
           min_eigv = 0) {
    stopifnot(inherits(y, "dgCMatrix"))
    irow = y@i
    pcol = y@p
    n <- dim(y)
    m <- n[2]
    n <- n[1]
    lambda.M <- lambda.M * m * n # rescalling to save time 
    k <- ncol(X)
    if (trace.it)
      nz = nnzero(y, na.counted = TRUE)
    #-------------------------------
    # the following part is for Laplacian regularization, only if similarity matrices are provided
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
    #' for warmup we extract  U, V, Dsq, and beta (named: u,d,v,beta)
    #' For initialization: Apply Naive MC followed by unregularized linear regression.
    warm = FALSE
    clean.warm.start(warm.start)
    if (!is.null(warm.start)) {
      #must have u,d and v components
      if (!all(match(c("u", "d", "v", "beta"), names(warm.start), 0) >
               0))
        stop("warm.start is missing some or all of the components.")
      warm = TRUE
      beta = warm.start$beta
      warm.out <- utils$prepare.M.warm.start(warm.start, J, n, m, min_eigv)
      U = warm.out$U
      V = warm.out$V
      Dsq = warm.out$Dsq
      
    } else{
      # initialize. Warm start is not provided
      Xterms = utils$GetXterms(X)
      Y_naive = naive_MC(as.matrix(y))
      beta = Xterms$X1 %*% Y_naive
      Xbeta <- X %*% beta   #svdH$u %*% (svdH$v  %*% Y_naive)
      M <- Y_naive - Xbeta
      M <- utils$svdopt(as.matrix(M), J, n, m, F, F)
      U = M$u
      V = M$v
      Dsq = M$d
      Y_naive <- Xbeta <- M <- NULL
      #---------------------------------------------------------------
      
    }
    
    #----------------------------------------
    yobs <- y@x # y will hold the model residuals
    ratio <- 1
    iter <- 0
    #----------------------------------------
    while ((ratio > thresh) & (iter < maxit)) {
      iter <- iter + 1
      U.old = U
      V.old = V
      Dsq.old = Dsq
      #----------------------------------------------
      # part 0: Update y (training residulas)
      # prereq: U, Dsq, V, Q, R, yobs
      # updates: y; VDsq; XQ
      VDsq = t(Dsq * t(V))
      M_obs = suvC(U, VDsq, irow, pcol)
      if (iter == 1) {
        xbeta.obs <- suvC(X, (t(beta)), irow, pcol)
        y@x = yobs - M_obs - xbeta.obs
      }
      #----------------------------------------------
      # part 1: update beta
      beta.old <- matrix(0, k, m)
      
      beta.thresh =   lambda.beta * learning.rate
      beta.iter = 0
      partial.update = learning.rate * as.matrix(t(X) %*% y + XtX %*% beta)
      while (sqrt(sum((beta - beta.old) ^ 2)) > 1e-3 &
             beta.iter < beta.iter.max) {
        beta.old <- beta
        update = beta + partial.update - learning.rate * XtX %*% beta
        beta <- matrix(0, k, m)
        beta[update > beta.thresh] = update[update > beta.thresh] - beta.thresh
        beta[update < -beta.thresh] = update[update < -beta.thresh] + beta.thresh
        beta.iter <- beta.iter + 1
        #print(sqrt(sum((beta - beta.old) ^ 2)))
      }
      
      #-------------------------------------------------------------------
      # part extra: re-update y
      xbeta.obs <-
        suvC(X, (t(beta)), irow, pcol)
      y@x = yobs - M_obs - xbeta.obs
      #--------------------------------------------------------------------
      # print("hi2")
      ##--------------------------------------------
      # part 3: Update B
      # prereq: U, VDsq, y, Dsq, lambda.M, L.b
      # updates U, Dsq, V
      B = crossprod(y, U) + VDsq # equiv to below
      # B = t(U) %*% y + t(VDsq)
      if (laplace.b)
        B = B - t(V)  %*% L.b
      D.star <- Dsq / (Dsq + lambda.M)
      B <- as.matrix(B %*% diag(D.star))
      #B = as.matrix(t((B) * (Dsq / (Dsq + lambda.M))))
      Bsvd = utils$svd_small_nc(B, FALSE, p = J) 
      V = Bsvd$u
      Dsq = Bsvd$d #pmax(Bsvd$d, min_eigv)
      U = U %*% (Bsvd$v)
      #-------------------------------------------------------------
      # part 4: Update A
      # prereq: U, D, VDsq, y, lambda.M, L.a
      # updates U, Dsq, V
      A = (y %*% V) + sweep(U, 2L, Dsq, `*`)
      #A = (y %*% V) + t(Dsq * t(U))
      if (laplace.a)
        A = A - L.a %*% U
      D.star <- Dsq / (Dsq + lambda.M)
      A <- sweep(A, 2L, D.star, `*`)
      #A = as.matrix(t(t(A) * (1 + lambda.M* Dsq^(-2))^2 ))
      Asvd =  utils$svd_small_nc(A, FALSE, p = J)
      U = Asvd$u
      Dsq = Asvd$d #pmax(Asvd$d, min_eigv)
      V = V %*% (Asvd$v)
      #------------------------------------------------------------------------------
      ratio =  utils$Frob(U.old, Dsq.old, V.old, U, Dsq, V)
      #------------------------------------------------------------------------------
      if (trace.it) {
        obj = (.5 * sum(y@x ^ 2) +
                 lambda.M * sum(Dsq)) / nz
        cat(iter,
            ":",
            "obj",
            format(round(obj, 5)),
            "ratio",
            ratio,
            "qiter",
            beta.iter,
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
    J = min(max(1,sum(Dsq > min_eigv)), J)
    
    out = list(
      u = U[, seq(J), drop = FALSE],
      d = Dsq[seq(J)],
      v = V[, seq(J), drop = FALSE],
      beta = beta,
      #---------------------
      # hyperparameters
      lambda.M = lambda.M,
      lambda.beta = lambda.beta,
      J = J,
      lambda.a = lambda.a,
      lambda.b = lambda.b,
      #--------------------------
      # convergence
      n_iter = iter
    )
    out
  }
