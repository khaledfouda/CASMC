# here we compare/test our implementation
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

# 0. prepare the data
print_performance <-
 function(inp,
          outp,
          test_error = error_metric$rmse,
          mao = FALSE,
          name = "",
          showres=T,
          rdig = 4) {
  error_function_name <-
   names(which(sapply(error_metric, identical, test_error)))
  
  if (mao) {
   M = outp$M
   O = outp$estimates
  } else{
   M = outp$u %*% (outp$d * t(outp$v))
   O = M + inp$X %*% outp$beta
  }
  
  error_beta <- error_M <- rank_beta <- rank_M <- NA
  
  if(!is.null(outp$beta))
    error_beta <- round(test_error(outp$beta, inp$beta), rdig)
  if(!is.null(M))
    error_M <- round(test_error(M, inp$M), rdig)
  error_test <-
   round(test_error(O[inp$W == 0], inp$O[inp$W == 0]), rdig)
  error_train <-
   round(test_error(O[inp$W != 0], inp$O[inp$W != 0]), rdig)
  error_Y <- 
   round(test_error(O[inp$W != 0], inp$Y[inp$W != 0]), rdig)
  
  if(!is.null(outp$beta))
    rank_beta <- qr(outp$beta)$rank
  if(!is.null(M))
    rank_M <- qr(M)$rank
  
  result_df <- data.frame(
   Metric = c(paste0(error_function_name, "(rank)")),
   model = name,
   Beta = sprintf(paste0("%.",rdig,"f(%2d)"), error_beta, rank_beta),
   M = sprintf(paste0("%.",rdig,"f(%2d)"), error_M, rank_M),
   test = sprintf(paste0("%.",rdig,"f"), error_test),
   train = sprintf(paste0("%.",rdig,"f"), error_train),
   Y = sprintf(paste0("%.",rdig,"f"), error_Y)
  )
  
  if(showres){
  print(knitr::kable(result_df, format = "simple" ))
  }else
 return(result_df)
 }


dat <-
 generate_simulation_rows(
  800,
  900,
  r = 10,
  k = 10, 
  missing_prob = 0.9,
  coll = F,
  prepare_for_fitting = TRUE,
  half_discrete = FALSE,
  informative_cov_prop = 0.7,
  mv_beta = T,
  seed = 2023
 )
X_r = reduced_hat_decomp(dat$X, 1e-2)


(900*10/800) * 10
svd(dat$X %*% naive_fit(as.matrix(dat$Y), dat$X)$beta)$d[1]


#' # 1
#+ something

fit_rank <- CASMC_cv_rank(
 y_train = dat$fit_data$train,
 X = X_r$X,
 y_valid = dat$fit_data$valid,
 W_valid = dat$fit_data$W_valid,
 y = dat$fit_data$Y,
 error_function = error_metric$rmse,
 lambda.factor = 1 / 4,
 lambda.init = NULL,
 n.lambda = 20,
 rank.init = 2,
 rank.limit = 30,
 rank.step = 2,
 pct = 0.98,
 lambda.a = 0,
 S.a = NULL,
 lambda.b = 0,
 S.b = NULL,
 early.stopping = 1,
 thresh = 1e-6,
 maxit = 30,
 trace = F,
 print.best = TRUE,
 quiet = FALSE,
 warm = NULL,
 rank_x = X_r$rank,
 r_min = 0,
 r_max = X_r$rank,
 track = TRUE,
 max_cores = 30,
 seed = 2023
)
fit_rank$fit$beta[,1:4]
dat$beta[,1:4]
fit_rank$rank_M
print_performance(dat, fit_rank$fit, error_metric$rmse, F, "CASMC(Rank)",F,3)


#' # 2
fit_l2 <- CASMC_cv_L2(
 y_train = dat$fit_data$train,
 X = X_r$X,
 y_valid = dat$fit_data$valid,
 W_valid = dat$fit_data$W_valid,
 y = dat$fit_data$Y,
 error_function = error_metric$rmse,
 lambda.factor = 1 / 4,
 lambda.init = NULL,
 n.lambda = 20,
 rank.init = 2,
 rank.limit = 30,
 rank.step = 2,
 pct = 0.98,
 lambda.a = 0,
 S.a = NULL,
 lambda.b = 0,
 S.b = NULL,
 early.stopping = 1,
 thresh = 1e-6,
 maxit = 100,
 trace = FALSE,
 print.best = TRUE,
 quiet = FALSE,
 warm = NULL,
 lambda.beta.grid = "default",
 track = TRUE,
 max_cores = 23,
 seed = 2023
)
fit_l2$fit$beta[,1:4]
dat$beta[,1:4]

print_performance(dat, fit_l2$fit, error_metric$rmse, F, "CASMC(L2)",F,3)

#' # 3
fit_l2_ <- CASMC_cv_L2(
 y_train = dat$fit_data$train,
 X = dat$X,
 y_valid = dat$fit_data$valid,
 W_valid = dat$fit_data$W_valid,
 y = dat$fit_data$Y,
 error_function = error_metric$rmse,
 lambda.factor = 1 / 4,
 lambda.init = NULL,
 n.lambda = 20,
 rank.init = 2,
 rank.limit = 30,
 rank.step = 2,
 pct = 0.98,
 lambda.a = 0,
 S.a = NULL,
 lambda.b = 0,
 S.b = NULL,
 early.stopping = 1,
 thresh = 1e-6,
 maxit = 2,
 trace = FALSE,
 print.best = TRUE,
 quiet = FALSE,
 warm = NULL,
 lambda.beta.grid = "default",
 track = TRUE,
 max_cores = 23,
 seed = 2023
)
#' # 4
simpute_fit <- simpute.cv(
  Y_train = as.matrix(dat$fit_data$train),
  y_valid = dat$fit_data$valid,
  W_valid = dat$fit_data$W_valid,
  y = dat$Y,
  n.lambda = 20,
  trace = FALSE,
  print.best = TRUE,
  tol = 5,
  thresh = 1e-6,
  rank.init = 2,
  rank.limit = 30,
  rank.step = 2,
  maxit = 300,
  seed= 2023

)

#' # 5
mao.fitB <- Mao.cv(
 Y = dat$Y,
 X = dat$X,
 W = dat$W,
 n_folds = 5,
 lambda.1_grid = seq(0, 1, length = 20),
 lambda.2_grid = seq(0.9, 0.1, length = 20),
 alpha_grid = c(1),#seq(0.992, 1, length = 5),
 seed = 2023,
 numCores = 1,
 n1n2_optimized = TRUE,
 test_error = error_metric$rmse,
 theta_estimator = Mao_weights$binomial,
 sequential = FALSE
)


#' # 6

mao.fit <- Mao.cv(
 Y = dat$Y,
 X = dat$X,
 W = dat$W,
 n_folds = 5,
 lambda.1_grid = seq(0, 1, length = 20),
 lambda.2_grid = seq(0.9, 0.1, length = 20),
 alpha_grid = c(1),#seq(0.992, 1, length = 5),
 seed = 2023,
 numCores = 1,
 n1n2_optimized = TRUE,
 test_error = error_metric$rmse,
 theta_estimator = Mao_weights$uniform,
 sequential = FALSE
)
#' # 7
mao.fitB$best_parameters

naive_model <- naive_fit(dat$Y, dat$X)


metric = error_metric$rmse
results <- rbind(
 print_performance(dat, fit_rank$fit, metric, F, "CASMC(Rank)",F,3),
 print_performance(dat, fit_l2$fit, metric, F, "CASMC(L2)",F,3),
 print_performance(dat, mao.fit$fit, metric, TRUE, "Mao(Uni)",F,3),
 print_performance(dat, mao.fitB$fit, metric, TRUE, "Mao(Binom)",F,3),
 print_performance(dat, simpute_fit, metric, TRUE, "SoftImpute", F, 3),
 #print_performance(dat, fit_rank_$fit, metric, F, "CASMC(Rank) 1 iter",F,3),
 print_performance(dat, fit_l2_$fit, metric, F, "CASMC(L2) 1 iter",F, 3),
 print_performance(dat, naive_model, metric, TRUE, "Naive", F, 3))

results  |> 
 arrange(Beta) |> 
 mutate(Beta = paste0("[",1:nrow(results),"]",Beta)) |> 
 arrange(M) |> 
 mutate(M = paste0("[",1:nrow(results),"]",M)) |> 
 arrange(Y) |> 
 mutate(Y = paste0("[",1:nrow(results),"]",Y)) |> 
 dplyr::select(-Y) |> 
 arrange(train) |> 
 mutate(train = paste0("[",1:nrow(results),"]",train)) |> 
 arrange(test) |> 
 mutate(test = paste0("[",1:nrow(results),"]",test))  ->
 results


kable( results, format = "simple")
#----------------------------------------------------------------------------------
#fit_l2_backup <- fit_l2
#fit_l2 <- fit_rank

fit_rank$fit$M <- unsvd(fit_rank$fit)
dat$xbeta = dat$X %*% dat$beta

fit_l2$fit$M <- unsvd(fit_l2$fit)
fit_l2$xbeta <- dat$X %*% fit_l2$fit$beta
fit_l2$O <- fit_l2$xbeta + fit_l2$fit$M

mao.fitB$xbeta <- dat$X %*%  mao.fitB$fit$beta
mao.fitB$O <- mao.fitB$xbeta +  mao.fitB$fit$M

abs(dat$xbeta[1,1:6] - mao.fitB$xbeta[1,1:6]) < abs(dat$xbeta[1,1:6] - fit_l2$xbeta[1,1:6])
abs(dat$M[1,1:6] - mao.fitB$fit$M[1,1:6]) < abs(dat$M[1,1:6] - fit_l2$fit$M[1,1:6])
abs(dat$O[1,1:6] - mao.fitB$O[1,1:6]) < abs(dat$O[1,1:6] - fit_l2$O[1,1:6])

#' #' # 8
#'
library(knitr)

comparison_df <- data.frame(
  Observation = 1:ncol(dat$xbeta),
  True_xbeta = dat$xbeta[1, ],
  Mao_FitB_xbeta = mao.fitB$xbeta[1, ],
  Fit_L2_xbeta = fit_l2$xbeta[1, ],
  True_M = dat$M[1, ],
  Mao_FitB_M = mao.fitB$fit$M[1, ],
  Fit_L2_M = fit_l2$fit$M[1, ],
  True_O = dat$O[1, ],
  Mao_FitB_O = mao.fitB$O[1, ],
  Fit_L2_O = fit_l2$O[1, ]
)

comparison_table <- comparison_df[1:10, c("Observation", "True_xbeta", "Mao_FitB_xbeta", "Fit_L2_xbeta",
                                      "True_M", "Mao_FitB_M", "Fit_L2_M",
                                      "True_O", "Mao_FitB_O", "Fit_L2_O")]
kable(comparison_table, format = "html", col.names = c("Obs", "True xbeta", "Mao FitB xbeta", "Fit L2 xbeta",
                                                       "True M", "Mao FitB M", "Fit L2 M",
                                                       "True O", "Mao FitB O", "Fit L2 O")) %>%
  kable_styling() %>%
  column_spec(2, background = "lightgreen") %>%
  column_spec(3, background = "lightblue") %>%
  column_spec(4, background = "lightblue") %>%
  column_spec(5, background = "lightgreen") %>%
  column_spec(6, background = "lightyellow") %>%
  column_spec(7, background = "lightyellow") %>%
  column_spec(8, background = "lightgreen") %>%
  column_spec(9, background = "lightpink") %>%
  column_spec(10, background = "lightpink")


#' # 9
#
#
library(ggplot2)
library(gridExtra)
ncols <- 100

# Create long format data for ggplot
long_comparison_df <- rbind(
  data.frame(Observation = 1:ncols, Value = dat$xbeta[1,1:100 ], Type = "True xbeta"),
  data.frame(Observation = 1:ncols, Value = mao.fitB$xbeta[1,1:100 ], Type = "Mao FitB xbeta"),
  data.frame(Observation = 1:ncols, Value = fit_l2$xbeta[1, 1:100], Type = "Fit L2 xbeta"),
  data.frame(Observation = 1:ncols, Value = dat$M[1, 1:100], Type = "True M"),
  data.frame(Observation = 1:ncols, Value = mao.fitB$fit$M[1,1:100 ], Type = "Mao FitB M"),
  data.frame(Observation = 1:ncols, Value = fit_l2$fit$M[1,1:100 ], Type = "Fit L2 M"),
  data.frame(Observation = 1:ncols, Value = dat$O[1,1:100 ], Type = "True O"),
  data.frame(Observation = 1:ncols, Value = mao.fitB$O[1,1:100 ], Type = "Mao FitB O"),
  data.frame(Observation = 1:ncols, Value = fit_l2$O[1,1:100 ], Type = "Fit L2 O")
)

# Create individual plots
plot_xbeta <- ggplot(long_comparison_df[long_comparison_df$Type %in% c("True xbeta", "Mao FitB xbeta", "Fit L2 xbeta"), ], aes(x = Observation, y = Value, color = Type, linetype = Type)) +
  geom_line() +
  scale_color_manual(values = c("True xbeta" = "black", "Mao FitB xbeta" = "blue", "Fit L2 xbeta" = "red")) +
  scale_linetype_manual(values = c("True xbeta" = "longdash", "Mao FitB xbeta" = "solid", "Fit L2 xbeta" = "solid")) +
  labs(title = "Comparison of xbeta", y = "Value", x = "Observation") +
  theme_minimal()

plot_M <- ggplot(long_comparison_df[long_comparison_df$Type %in% c("True M", "Mao FitB M", "Fit L2 M"), ], aes(x = Observation, y = Value, color = Type, linetype = Type)) +
  geom_line() +
  scale_color_manual(values = c("True M" = "black", "Mao FitB M" = "blue", "Fit L2 M" = "red")) +
  scale_linetype_manual(values = c("True M" = "longdash", "Mao FitB M" = "solid", "Fit L2 M" = "solid")) +
  labs(title = "Comparison of M", y = "Value", x = "Observation") +
  theme_minimal()

plot_O <- ggplot(long_comparison_df[long_comparison_df$Type %in% c("True O", "Mao FitB O", "Fit L2 O"), ], aes(x = Observation, y = Value, color = Type, linetype = Type)) +
  geom_line() +
  scale_color_manual(values = c("True O" = "black", "Mao FitB O" = "blue", "Fit L2 O" = "red")) +
  scale_linetype_manual(values = c("True O" = "longdash", "Mao FitB O" = "solid", "Fit L2 O" = "solid")) +
  labs(title = "Comparison of O", y = "Value", x = "Observation") +
  theme_minimal()

# Combine the plots into one graph with 3 rows and 1 column
grid.arrange(plot_xbeta, plot_M, plot_O, nrow = 3)
# 


set.seed(42)  # For reproducibility

# User and Movie Data
nUsers <- 30
nMovies <- 40
ageGroups <- c("young", "middle_aged", "senior")
genres <- c("action", "comedy", "drama", "sci-fi")

Users <- matrix(sample(ageGroups, size = nUsers, replace = TRUE), ncol = 1)
Movies <- matrix(sample(genres, size = nMovies, replace = TRUE), ncol = 1)

# Ratings
Yraw <- matrix(0, nrow = nUsers, ncol = nMovies)
for (i in 1:nUsers) {
  moviesToRate <- sample(nMovies, size = 20, replace = FALSE)
  Yraw[i, moviesToRate] <- sample(1:5, size = 20, replace = TRUE)
}
R <- Matrix::Matrix(Yraw != 0, sparse = TRUE)

# Features
featureIDs <- new.env()  # Use environment for dictionary-like functionality
currentID <- 0
for (age in ageGroups) {
  for (genre in genres) {
    featureIDs[[paste(age, genre)]] <- currentID
    currentID <- currentID + 1
  }
}

flatX <- matrix(0, nrow = nUsers, ncol = nMovies)
for (i in 1:nUsers) {
  for (j in 1:nMovies) {
    flatX[i, j] <- featureIDs[[paste(Users[i, ], Movies[j, ])]]
  }
}

# Proportions
Ps <- sapply(featureIDs, function(id) {
  mean(flatX == id)
})

# Convert to Torch tensors if needed (using the 'torch' package)
# Ps <- torch::torch_tensor(Ps)
# flatX <- torch::torch_tensor(flatX)


# flatX is a 2d Matrix

categs <- sort(unique(as.vector(flatX)))

Ps <- vector(length=length(categs))
for(i in 1:length(categs))
  Ps[i] <- mean(flatX==categs[i])


