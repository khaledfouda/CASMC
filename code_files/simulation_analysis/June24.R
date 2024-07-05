dat <-
 generate_simulation_rows(
  800,
  700,
  r = 10,
  k = 20, 
  missing_prob = 0.9,
  coll = F,
  prepare_for_fitting = TRUE,
  half_discrete = FALSE,
  informative_cov_prop = .1,mar_sparse = F,
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
 M_cv_param = list(
   rank.init = 2,
   rank.limit = 30,
   rank.step = 2,
   pct = 0.98,
   lambda.factor = 1/4,
   lambda.init = NULL,
   n.lambda = 20, 
   early.stopping = 1
 ),
 beta_cv_param = list(
   rank.init = 2,
   rank.limit = qr(dat$X)$rank,
   rank.step = 2,
   pct = 0.98,
   lambda.multi.factor = 20,
   lambda.init = NULL,
   n.lambda = 20, 
   early.stopping = 1
 ),
 use_warmstart = TRUE,
 quiet = F,
 trace = F,  
 track = F, step3 = T,
 seed = 2023
)

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
##########################################################
residuals2 <- dat$O[dat$W==0] -  fit.$estimates[dat$W==0]


# 1. residual analysis
residual_analysis(residuals2)

# 2. checking the likelihood
logLikelihood(residuals)
Likelihood_ratio_index(logLikelihood(residuals), logLikelihood(residuals2))
Cox_Snell_R2(logLikelihood(residuals), logLikelihood(residuals2), length(residuals))


LogLik0 <- logLikelihood(residuals2)

SImpute_Sim_Wrapper(dat)
Mao_Sim_Wrapper(dat, LogLik_SI = LogLik0)
CASMC_0_Sim_Wrapper(dat, LogLik_SI = LogLik0)

prepare_output(Sys.time(), fit.$estimates, dat$O, dat$W)
###############################################################
start_time = Sys.time() 
fitkf <- CASMC2_cv_kf(
  Y = dat$Y,
  X = dat$X,
  obs_mask = dat$W,
  y = dat$fit_data$Y,
  error_function = error_metric$rmse,
  warm = NULL,
  M_cv_param = list(
    rank.init = 2,
    rank.limit = 30,
    rank.step = 2,
    pct = 0.98,
    lambda.factor = 1/4,
    lambda.init = NULL,
    n.lambda = 20, 
    early.stopping = 1
  ),
  beta_cv_param = list(
    rank.init = 2,
    rank.limit = qr(dat$X)$rank,
    rank.step = 2,
    pct = 0.98,
    lambda.multi.factor = 20,
    lambda.init = NULL,
    n.lambda = 20, 
    early.stopping = 1
  ),
  quiet = F,
  trace = F,  
  track = F,
  max_cores = 20,
  n_folds = 5,
  seed = 2023,
)

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
set.seed(2023); CASMC_3a_Sim_Wrapper(dat, max_cores = 20) 


system.time({
  for (i in 1:5000)
    x = 1:10 * diag(1,10,10)
  
})

x = kronecker(diag(1:10), matrix(5,10,10))
x[1:20,1:20]



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








