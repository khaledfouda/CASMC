dat <-
 generate_simulation_rows(
  800,
  700,
  r = 10,
  k = 13, 
  missing_prob = 0.9,
  coll = F,
  prepare_for_fitting = TRUE,
  half_discrete = FALSE,
  informative_cov_prop = .5,mar_sparse = F,
  mv_beta = T,
  seed = 2023
 )



start_time = Sys.time()
fiti <- CASMC2_cv2(
 y_train = dat$fit_data$train,
 X = dat$X,
 y_valid = dat$fit_data$valid,
 W_valid = dat$fit_data$W_valid,
 y = dat$fit_data$Y,
 error_function = error_metric$rmse,
 warm = NULL,
 quiet = F,
 trace = T,  
 track = T, step3 = T,
 #rank.beta.init = 10, rank.beta.limit = 10, lambda.beta.grid = c(0,0),
 rank.beta.step = 1, early.stopping = 10,
 lambda.beta.length = 80,
 lambda.beta.grid = "default1",
 max_cores = max_cores,
 seed = 2023,
)



# set.seed(2023); fit. <- CASMC2_fit(
#   y = dat$fit_data$Y,
#   X = dat$X,
#   J = fiti$hparams$rank.M,
#   lambda.M = fiti$hparams$lambda.M,
#   r = fiti$hparams$rank.beta,
#   lambda.beta = fiti$hparams$lambda.beta,
#   warm.start = NULL,
#   trace.it = F, 
#   thresh = 1e-6,
#   maxit = 300
# )
# 
# test_error <- error_metric$rmse_normalized
# fit.$beta <- list(u=fit.$ub, d=fit.$db^2, v=fit.$vb)

fit. = fiti$fit
fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
fit.$beta = unsvd(fit.$beta)
fit.$estimates = fit.$M + dat$X %*% fit.$beta

results = list(model = "CASMC-2")
results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
results$lambda.M = fiti$hparams$lambda.M
results$lambda.beta = fiti$hparams$lambda.beta
results$error.test = test_error(fit.$estimates[dat$W == 0], dat$O[dat$W == 0])
results$error.all = test_error(fit.$estimates, dat$O)
results$error.M = test_error(fit.$M, dat$M)
results$error.beta = test_error(fit.$beta, dat$beta)
results$rank_M = sum(fit.$d > 0)
results$rank_beta = qr(fit.$beta)$rank
results$sparse_in_sparse = sum(dat$beta == 0 & fit.$beta == 0) /
  (sum(dat$beta == 0) +  1e-17)
results$sparse_in_nonsparse = sum(dat$beta != 0 &
                                    fit.$beta == 0) /
  (sum(dat$beta != 0) +  1e-17)
results
####################################################
start_time = Sys.time() 
fitkf <- CASMC2_cv_kf(
  Y = dat$Y,
  X = dat$X,
  obs_mask = dat$W,
  y = dat$fit_data$Y,
  error_function = error_metric$rmse,
  warm = NULL,
  quiet = F,
  trace = F,  
  track = F,
  #rank.beta.init = 10, rank.beta.limit = 10, lambda.beta.grid = c(0,0),
  rank.beta.step = 1, early.stopping = 10,
  lambda.beta.length = 80,
  lambda.beta.grid = "default1",
  max_cores = 20,
  n_folds = 10,
  seed = 2023,
)

test_error <- error_metric$rmse_normalized

fitkf$estimates = fitkf$M + dat$X %*% fitkf$beta

results = list(model = "CASMC-2")
results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
results$lambda.M = fitkf$M_param$lambda
results$lambda.beta = fitkf$beta_param$lambda
results$error.test = test_error(fitkf$estimates[dat$W == 0], dat$O[dat$W == 0])
results$error.all = test_error(fitkf$estimates, dat$O)
results$error.M = test_error(fitkf$M, dat$M)
results$error.beta = test_error(fitkf$beta, dat$beta)
results$rank_M = qr(fitkf$M)$rank
results$rank_beta = qr(fitkf$beta)$rank
results$sparse_in_sparse = sum(dat$beta == 0 & fitkf$beta == 0) /
  (sum(dat$beta == 0) +  1e-17)
results$sparse_in_nonsparse = sum(dat$beta != 0 &
                                    fitkf$beta == 0) /
  (sum(dat$beta != 0) +  1e-17)
results

####################################################
set.seed(2023); CASMC_0_Sim_Wrapper(dat, max_cores = 20) 



fit.$beta[1:5,1:5]
dat$beta[1:5,1:5]


system.time({for(i in 1:2) fit. = CASMC2_fit(
  y = dat$fit_data$train,
  X = dat$X,
  J = 4,
  lambda.M = 2,
  r = 4,
  lambda.beta = 1,
  warm.start = NULL,
  trace.it = F
)})


#---------------------------------------------------------------

k = 10
r = 6
Q =  matrix(rnorm(k*r), k, r)
all(Q==matrix(as.vector(Q), k, r))

kk = kronecker(diag(c(1:r),r,r),matrix(rnorm(k*k),k,k))  
kk[kk==0] <- NA
kk <- as(kk, "Incomplete")

length(kk@x)

all( kk %*% as.vector(Q) == kk %*% matrix(as.vector(Q),ncol = 1))


as.vector(Q) + kk %*% as.vector(Q)

all(as.vector(3 * Q) ==  diag(3, 60,60) %*% as.vector(Q)) 

#-----------------------------------------------------------
# 1. fixed part of Q as a vector -> P1V
# 2. get kronecker(XtX, diag(Dsq)) ; save as Incomplete -> XtXD
# 3. get Hessian as: XtXD + diag(lambda, kr, kr)
# 4. set  as.vector(Q) -> Q
# 5. Start Loop:
# 6. Gradient(as vec) = P1V + (XtXD + diag(lambda, kr, kr)) Q 
# 7. Update:  Q <- Q - solve(H(Q)) %*% Grad(Q)
#-----------------------------------------------------------








