setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
# source("./code_files/import_lib.R")

load_model_bixi_dat2 <- function(time_cov=TRUE, seed = 2023, testp=0.2,
                                 missp = 0.2, validp = 0.2){

  if(!is.null(seed)) set.seed(seed)

bdat <- BixiData$new()
i = ifelse(time_cov, 2, 1)

model.dat <- list()

train.df <-
  readRDS(paste0("./BIXI/data/splits/split_", i, "_train.rds"))
test.df <-
  readRDS(paste0("./BIXI/data/splits/split_", i, "_test.rds"))

if(time_cov){
train.df |>
  as.data.frame() |>
  mutate(time = as.Date(time)) %>%
  arrange(location, time) ->  train.df

test.df |>
  as.data.frame() |>
  mutate(time = as.Date(time)) %>%
  arrange(location, time) ->  test.df

  
}else{
  train.df |>
    as.data.frame() |>
    mutate(time = as.Date(time)) %>%
    rename(location = time, time = location) %>%
    arrange(location, time) ->  train.df
  
  test.df |>
    as.data.frame() |>
    mutate(time = as.Date(time)) %>%
    rename(location = time, time = location) %>%
    arrange(location, time) ->  test.df
  
}
dim(train.df)
dim(test.df)
combined <- 
rbind(train.df, test.df) %>%
  group_by(location, time) %>% 
  filter(!(is.na(nb_departure) & any(!is.na(nb_departure)))) %>% 
  ungroup() ->
  combined
dim(combined) 

sum(is.na(train.df$nb_departure)) / nrow(train.df)
sum(is.na(combined$nb_departure)) / nrow(combined)
nrow(test.df)



X <- combined |>
  group_by(time) |>
  filter(row_number() == 1) |>
  ungroup() |>
  select(-location, -time, -nb_departure)

dim(X)
model.dat$X <- as.matrix(X)

reshape2::dcast(combined, time ~ location, value.var = "nb_departure") |>
  select(-time) |>
  as.matrix() ->
  Y
dim(Y)
model.dat$Y <- Y
sum(is.na(as.vector(Y)))


obs.mask <- (1 * !is.na(Y)) %>%
  {
    colnames(.) <- NULL
    .
  } %>%
  as.matrix()
print(sum(obs.mask == 1) / length(obs.mask))
#--------------------------------------------------
cor_matrix <- cor(naive_MC(Y), X)
tau <- 0.4
high_corr_indices <- which(apply(abs(cor_matrix), 1, max) > tau)
# high_corr_indices <- c(
#                         sample(1:nrow(Y),
#                                length(high_corr_indices), replace=FALSE))
m <- nrow(Y)
missing_rate <- .95
test_mask <- matrix(1, nrow(Y), ncol(Y))
for (j in high_corr_indices) {
  missing_indices <- sample(1:m, floor(missing_rate * m), replace = FALSE)
  test_mask[missing_indices, j] <- 0
}
test_mask[obs.mask==0] = 1
sum(test_mask==0) / length(Y)
#--------------------------------------------------
# divide into train and test.
masks <- list(obs = obs.mask)
masks$test = test_mask
#masks$test = utils$MC_train_test_split(obs.mask, testp, seed)
obs = masks$obs * masks$test
masks$miss = utils$MC_train_test_split(obs, missp, seed)
obs = obs * masks$miss

#masks$test = masks$test * masks$miss
masks$valid = utils$MC_train_test_split(obs, validp, seed)

sum(masks$test==0) / length(masks$test)
sum(masks$valid==0) / length(masks$test)
sum(masks$valid==0 & masks$obs==0) / length(masks$test)

fit_data <- list()
fit_data$W_valid <- masks$valid
fit_data$train <- utils$to_incomplete(Y * masks$valid * masks$test * masks$miss)
fit_data$valid <- Y[masks$valid== 0]
fit_data$Y = utils$to_incomplete(Y * masks$test * masks$miss)

out <- list(
  O = Y * masks$miss,
  W = masks$test,
  X = as.matrix(X),
  Y = Y * masks$test* masks$miss,
  fit_data = fit_data,
  masks = masks
)
return(out)
# 
# train.df |>
#   select(time, location, nb_departure) |>
#   merge(
#     select(test.df, time, location, nb_departure),
#     by = c("time", "location"),
#     all.x = T
#   ) |>
#   as.data.frame() |>
#   arrange(location, time)  |>
#   mutate(missing =  !(is.na(nb_departure.x) &
#                         (!is.na(nb_departure.y)))) -> mixed
# 
# print(sum(1 - mixed$missing) / nrow(mixed))
# 
# reshape2::dcast(mixed, time ~ location, value.var = "missing") %>%
#   select(-time) %>%
#   {
#     colnames(.) <- NULL
#     .
#   } %>%
#   as.matrix() ->
#   #as.integer() ->
#   test.mask
# 
# test.mask = 1 * test.mask
# 
# print(sum(1 - test.mask) / length(test.mask))
# obs.mask[1:5, 1:5]
# test.mask[1:5, 1:5]
# 
# valid_mask <- utils$MC_train_test_split(obs.mask, testp = 0.2)
# 
# 
# model.dat$masks <-
#   list(tr_val = obs.mask,
#        test = test.mask,
#        valid = valid_mask)
# #-----------------------
# # get test matrix
# train.df |>
#   select(time, location, nb_departure) |>
#   merge(
#     select(test.df, time, location, nb_departure),
#     by = c("time", "location"),
#     all.x = T
#   ) |>
#   as.data.frame() |>
#   arrange(location, time)  |>
#   select("time", "location", nb_departure.y) %>%
#   reshape2::dcast(time ~ location, value.var = "nb_departure.y") |>
#   select(-time) |>
#   as.matrix() %>%
#   utils$to_incomplete() ->
#   test
# 
# model.dat$splits <- list(
#   train = utils$to_incomplete(Y * valid_mask),
#   valid = utils$to_incomplete(Y * (1 - valid_mask)),
#   #Y[valid_mask==0],
#   test = test,
#   #Y = Y,
#   Y = utils$to_incomplete(Y)
#   
# )
# print(length(model.dat$splits$train@x) / length(Y))
# print(length(model.dat$splits$test@x))
# print(model.dat$splits$valid@x %>% length)
#-------------------------------
# get all observed matrix
# train.df |>
#   select(time, location, nb_departure) |>
#   merge(
#     select(test.df, time, location, nb_departure),
#     by = c("time", "location"),
#     all.x = T
#   ) |>
#   as.data.frame() |>
#   arrange(location, time)  |>
#   select("time", "location", nb_departure.x, nb_departure.y) %>%
#   mutate(nb_departure = ifelse(is.na(nb_departure.x), nb_departure.y, nb_departure.x)) %>%
#   reshape2::dcast(time ~ location, value.var = "nb_departure") |>
#   select(-time) |>
#   as.matrix() %>%
#   to_incomplete() ->
#   observed
# 
# print(length(observed@x))
# model.dat$depart <- Y#observed

# return(model.dat)
}
#--------------------------------------------------------------------------------------------
# length(model.dat$splits$Y@x)
# sum(model.dat$masks$tr_val!=0)
# length(model.dat$splits$train@x)
# 
# sum(!is.na(model.dat$Y))
# 
# model.dat <- dat <-  load_model_bixi_dat()
# SImpute_Bixi_Wrapper(model.dat)
# Mao_Bixi_Wrapper(model.dat)
# CASMC_0_Bixi_Wrapper(model.dat, train_on_all = TRUE)
# 
# CASMC_2_Bixi_Wrapper(model.dat, train_on_all = TRUE)
# 
# 
# CASMC_3a_Bixi_Wrapper(model.dat, train_on_all = TRUE) -> results
# 
# 
# Naive_Bixi_Wrapper(model.dat)
