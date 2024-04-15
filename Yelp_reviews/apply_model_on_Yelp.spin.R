
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")
source("./Yelp_reviews/load_yelp_data.R")
#---------------
print_rmse <- function(X_test, X_hat, model_name) {
 rmse <- sqrt(mean((X_test@x - X_hat@x) ^ 2))
 cat(sprintf("RMSE for %s is: %.4f\n", model_name, rmse))
}




scale = F
dat <- load_Yelp_data(scale=T, seed=2024)


best_fit = CASMC_cv_holdout_v2(
 dat$valid$train,
 dat$X_r,
 dat$valid$valid@x,
 dat$valid$W,
 dat$train.inc,
 trace = F,
 thresh = 1e-6,
 n.lambda = 30,
 rank.step = 2,
 rank.limit = 30 
)
fit1 = best_fit$fit
beta =  fit1$beta
M = fit1$u %*% (fit1$d * t(fit1$v))
A = M + dat$X_r$X %*% t(beta)
if(scale) A = revertBiScaledMatrix(as.matrix(A), dat$biScale)
#qr(A)$rank
preds <- A[dat$W_test==0]
RMSE_error(preds, dat$test.inc@x)

casmc_error = RMSE_error(dat$test.inc@x, preds)
#------------------------------------------------------------
dat$valid$train.mat <- as.matrix(dat$valid$train)
dat$train.mat <- as.matrix(dat$train.inc)
dat$valid$train.mat[dat$valid$train.mat==0] = NA
dat$train.mat[dat$train.mat==0] = NA
# soft-impute
start_time = Sys.time()
sout <- simpute.cv(
 dat$valid$train.mat,
 dat$train.mat,
 trace = FALSE,
 rank.limit = 30,
 print.best = FALSE,
 rank.step = 2
)
if(scale) sout$estimates <-  revertBiScaledMatrix(as.matrix(sout$estimates), dat$biScale)
preds_simpute <- sout$estimates[dat$W_test==0]
simpute_error = RMSE_error(dat$test.inc@x, preds_simpute)


time_softImpute = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
#---------------------------------------------------------------------------------
# Naive
start_time = Sys.time()
estimates = naive_MC(dat$train.mat)
if(scale) estimates =  revertBiScaledMatrix(as.matrix(estimates), dat$biScale)
preds_naive = expm1(estimates[dat$W_test==0])
time_naive <-
 as.numeric(difftime(Sys.time(), start_time, units = "secs"))
naive_error = RMSE_error(dat$test.inc@x, preds_naive)
naive_error
#----------------------------------------------------------------------------------




# source("./code_files/import_lib.R")

# dat <- load_Yelp_data(scale=FALSE, seed=2024)

#---------------------
library(cmfrec)
library(Matrix)
library(MatrixExtra)


X_train <- as.coo.matrix(dat$train.inc)
str(X_train)
X_test <- as.coo.matrix(dat$test.inc)
str(X_test)

model.classic <-
 CMF(
  X_train,
  k = 25,
  lambda = 0.1,
  scale_lam = TRUE,
  verbose = FALSE,
  nthreads = 6
 )



pred_classic <- predict(model.classic, X_test)
print_rmse(X_test, pred_classic, "classic model")

model.baseline <- MostPopular(X_train, lambda = 10, scale_lam = FALSE)
pred_baseline <- predict(model.baseline, X_test)
print_rmse(X_test, pred_baseline, "non-personalized model")

model.improved <- CMF(
 X_train,
 k = 25,
 lambda = 0.1,
 scale_lam = TRUE,
 add_implicit_features = TRUE,
 w_main = 0.75,
 w_implicit = 0.25,
 use_cg = FALSE,
 niter = 30,
 verbose = FALSE,
 nthreads = 6
)
pred_improved <- predict(model.improved, X_test)
print_rmse(X_test, pred_improved, "improved classic model")

model.w.sideinfo <- CMF(
 X_train,
 U = dat$X_r$X,
 k = 25,
 lambda = 0.1,
 scale_lam = TRUE,
 niter = 30,
 use_cg = FALSE,
 include_all_X = FALSE,
 w_main = 0.75,
 w_user = 0.5,
 w_item = 0.5,
 w_implicit = 0.5,
 center_U = FALSE,
 center_I = FALSE,
 nthreads = 6,
 verbose = FALSE
)
pred_side_info <- predict(model.w.sideinfo, X_test)
print_rmse(X_test, pred_side_info, "model with side info")

#------------------------------------------------
results <- data.frame(
 NonPersonalized = RMSE_error(X_test@x, pred_baseline@x),
 ClassicalModel = RMSE_error(X_test@x, pred_classic@x),
 ClassicPlusImplicit = RMSE_error(X_test@x, pred_improved@x),
 CollectiveModel = RMSE_error(X_test@x, pred_side_info@x),
 CASMC = casmc_error,
 SoftImpute = simpute_error,
 Naive = naive_error
)
results <- as.data.frame(t(results))
names(results) <- "RMSE"
results %>%
 kable() %>%
 kable_styling()

detach("package:cmfrec", unload = TRUE)
detach("package:MatrixExtra", unload = TRUE)
