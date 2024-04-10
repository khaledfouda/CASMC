
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
suppressMessages(suppressWarnings(source("./code_files/import_lib.R", local = FALSE)))
source("./MovieLens/load_data.R")
#---------------
print_rmse <- function(X_test, X_hat, model_name) {
 rmse <- sqrt(mean((X_test@x - X_hat@x) ^ 2))
 cat(sprintf("RMSE for %s is: %.4f\n", model_name, rmse))
}





dat <- load_movielens_100k("a", remove_bad_movies = T, scale = T)

best_fit = CASMC_cv_holdout_v2(
  dat$valid$train,
  dat$X_r,
  dat$valid$test,
  dat$valid$W,
  dat$train.mat,
  trace = F,
  thresh = 1e-6,
  n.lambda = 30,
  rank.step = 4,
  rank.limit = 50 
)
fit1 = best_fit$fit
beta =  fit1$beta
M = fit1$u %*% (fit1$d * t(fit1$v))
A = M + dat$X_r$X %*% t(beta)
A = revertBiScaledMatrix(as.matrix(A), dat$biScale)
#qr(A)$rank
preds <- A[cbind(dat$test.df$user_id, dat$test.df$movie_id)]
RMSE_error(preds, dat$test.df$rating)
#------------------------------------------------------------




#---------------------
library(cmfrec)
library(Matrix)
library(MatrixExtra)

dat <-
 load_movielens_100k("a", remove_bad_movies = TRUE, scale = FALSE)








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
 CASMAC = RMSE_error(dat$test.df$rating, preds)
)
results <- as.data.frame(t(results))
names(results) <- "RMSE"
results %>%
 kable() %>%
 kable_styling()
