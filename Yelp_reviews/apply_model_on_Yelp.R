
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
#---------------
print_rmse <- function(X_test, X_hat, model_name) {
 rmse <- sqrt(mean((X_test@x - X_hat@x) ^ 2))
 cat(sprintf("RMSE for %s is: %.4f\n", model_name, rmse))
}


all_results = list()
covariates <- "rows"
for(covariates in c("rows", "columns")){
source("./code_files/import_lib.R")
source("./Yelp_reviews/load_yelp_data.R")
scale = F
dat <- load_Yelp_data(scale=T, seed=2024,covariates = covariates,subset = "_4x3")


best_fit = CASMC_cv_holdout_with_r(
 dat$valid$train,
 dat$X_r,
 dat$valid$valid@x, 
 dat$valid$W,
 r_min = 0,
 y = dat$train.inc,
 trace = F,
 thresh = 1e-6,
 n.lambda = 30,
 rank.step = 2,
 rank.limit = 30,
 track_r = T,
 max_cores = ceiling( (dat$X_r$rank+2)/2)
)
fit1 = best_fit$fit
beta =  fit1$Beta$u %*% (fit1$Beta$d * t(fit1$Beta$v))
M = fit1$u %*% (fit1$d * t(fit1$v))
A = M + dat$X_r$X %*% t(beta)
if(scale) A = revertBiScaledMatrix(as.matrix(A), dat$biScale)
#qr(A)$rank
preds <- A[dat$W_test==0]
RMSE_error(preds, dat$test.inc@x)
print(best_fit$r)
casmc_error = RMSE_error(dat$test.inc@x, preds)
#------------------------------------------------------------
dat$valid$train.mat <- as.matrix(dat$valid$train)
dat$train.mat <- as.matrix(dat$train.inc)
dat$valid$train.mat[dat$valid$train.mat==0] = NA
dat$train.mat[dat$train.mat==0] = NA
# soft-impute
start_time = Sys.time()
# sout <- simpute.cv(
#  dat$valid$train.mat,
#  dat$train.mat,
#  trace = FALSE,
#  rank.limit = 30,
#  print.best = FALSE,
#  rank.step = 2
# )
# if(scale) sout$estimates <-  revertBiScaledMatrix(as.matrix(sout$estimates), dat$biScale)
# preds_simpute <- sout$estimates[dat$W_test==0]
simpute_error = 4#RMSE_error(dat$test.inc@x, preds_simpute)


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
  lambda = 0.1, #user_bias = F, item_bias = F,
  scale_lam = TRUE,
  verbose = FALSE,
  nthreads = 6
 )



pred_classic <- predict(model.classic, X_test)
print_rmse(X_test, pred_classic, "classic model")

model.baseline <- MostPopular(X_train, lambda = 10, scale_lam = FALSE)#,user_bias = T)
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
 niter = 30,# user_bias = F, item_bias = F,
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
 #user_bias = F, item_bias = F,
 verbose = FALSE
)
pred_side_info <- predict(model.w.sideinfo, X_test)
print_rmse(X_test, pred_side_info, "model with side info")

detach("package:cmfrec", unload = TRUE)
detach("package:MatrixExtra", unload = TRUE)
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


all_results[[covariates]] <- 
  results %>%
 kable(format = "pipe") %>%
 kable_styling() #%>%
  #print()

}


all_results[["rows"]]

all_results[["columns"]]
