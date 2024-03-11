# Soft Impute with Covariates fit. Fixed Lambda2. beta_partial is input.
simpute.svd.cov <-
   function (Y,
             X,
             beta_partial,
             J = 2,
             thresh = 1e-05,
             lambda = 0,
             maxit = 100,
             trace.it = FALSE,
             warm.start = NULL,
             min_rank = 2,
             svd_maxit = 100) {
      # are you scaling???
      
      n1 <- dim(Y)[1]
      n2 <- dim(Y)[2]
      m1 <- dim(X)[2]
      ynas <- Y == 0
      nz = n1 * n2 - sum(ynas)
      J = max(J, min_rank)
      yfill <- Y
      yfill[ynas] <- 0
      beta.estim <- beta_partial %*% yfill
      Xbeta <- X %*% beta.estim
      yplus <- yfill - Xbeta
      
      if (!is.null(warm.start)) {
         ###must have u,d and v components
         warm.start = clean.warm.start(warm.start)
         if (!all(match(c("u", "d", "v", "beta.estim"), names(warm.start), 0) >
                  0))
            stop("warm.start does not have components u, d and v")
         D = warm.start$d
         nzD = sum(D > 0)
         JD = min(nzD, J)
         U = warm.start$u[, seq(JD), drop = FALSE]
         V = warm.start$v[, seq(JD), drop = FALSE]
         D = D[seq(JD)]
         yhat = U %*% (D * t(V))
         Xbeta <- X %*% warm.start$beta.estim
         yfill[ynas] <- yhat[ynas] + Xbeta[ynas]
         beta.estim <- beta_partial %*% yfill
         Xbeta <- X %*% beta.estim
         yplus <- yfill - Xbeta
      }
      
      svd.yfill = propack.svd(yplus, J, opts = list(kmax = svd_maxit)) #svd(yplus) **
      ratio <- 1
      iter <- 0
      while ((ratio > thresh) & (iter < maxit)) {
         iter <- iter + 1
         svd.old = svd.yfill
         d = svd.yfill$d
         d = pmax(d - lambda, 0)
         # update
         yhat <-
            svd.yfill$u %*% (d * t(svd.yfill$v))
         # ** svd.yfill$u[, seq(J)] %*% (d[seq(J)] * t(svd.yfill$v[, seq(J)]))
         yfill[ynas] <- yhat[ynas] + Xbeta[ynas]
         # new svd
         beta.estim <- beta_partial %*% yfill
         Xbeta <- X %*% beta.estim
         yplus <- yfill - Xbeta
         svd.yfill = propack.svd(yplus, J, opts = list(kmax = svd_maxit)) #svd(yplus)
         
         #-- performance check
         ratio = Frob(
            svd.old$u[, seq(J)],
            d[seq(J)],
            svd.old$v[, seq(J)],
            svd.yfill$u[, seq(J)],
            pmax(svd.yfill$d - lambda, 0)[seq(J)],
            svd.yfill$v[, seq(J)]
         )
         if (trace.it) {
            obj = (.5 * sum((yfill - yhat)[!ynas] ^ 2) + lambda * sum(d)) / nz
            cat(iter, ":", "obj", format(round(obj, 5)), "ratio", ratio, "\n")
         }
      }
      # ** d = pmax(svd.yfill$d[seq(J)] - lambda, 0)
      d = pmax(svd.yfill$d - lambda, 0)
      # ** J = min(sum(d > 0) + 1, J)
      J = sum(d > 0) 
      svd.yfill = list(
         u = svd.yfill$u[, seq(J)],
         d = d[seq(J)],
         v = svd.yfill$v[, seq(J)],
         lambda = lambda,
         beta.estim = beta.estim
      )
      if (iter == maxit &
          trace.it)
         warning(paste("Convergence not achieved by", maxit, "iterations"))
      svd.yfill
   }
