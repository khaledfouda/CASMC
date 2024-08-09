setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
# source("./code_files/import_lib.R")

load_model_bixi_dat3 <- function(time_cov=TRUE, seed = 2023,
                                 validp = 0.2){

  if(!is.null(seed)) set.seed(seed)

#bdat <- BixiData$new()
i = 1 #ifelse(time_cov, 2, 1)

out <- list()

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


#-------------------------------------------------------
# needs:
#'  Y: matrix of observed
#'  X: Covariates matrix
#'  A: matrix of observed + test
#'  W: mask matrix with 0 for test
#'  model.dat:
#'    Y: sparse of Y
#'    train: sparse of Y * W_valid
#'    W_valid: mask with 0 for valid
#'    valid: vector of validation response.
#----------------------------------------------------------
train.df %>%
  arrange(location, time) %>% 
  group_by(location, time) %>% 
reshape2::dcast(time ~ location, value.var = "nb_departure") |>
  select(-time) |>
  as.matrix() ->
  Y
colnames(Y) <- NULL

out$Y <- Y

rbind(train.df, test.df) %>%
  group_by(location, time) %>% 
  filter(!(is.na(nb_departure) & any(!is.na(nb_departure)))) %>% 
  ungroup() %>% 
  arrange(location, time) %>% 
  reshape2::dcast(time ~ location, value.var = "nb_departure")  %>% 
  select(-time)  %>% 
  as.matrix() ->
  A 
colnames(A) <- NULL  
out$A <- out$O <- A
out$W <- !(is.na(Y) & (!is.na(A)))
print(sum(out$W==0) / length(out$W))
print(sum(is.na(Y)) / length(out$W)-(sum(is.na(A)) / length(out$W)) )
print(sum(is.na(A)) / length(out$W))
print(sum(is.na(Y)) / length(out$W))

#-- get X
rbind(train.df, test.df) %>%
  group_by(location, time) %>% 
  filter(!(is.na(nb_departure) & any(!is.na(nb_departure)))) %>% 
  ungroup() %>% 
  arrange(location, time) %>% 
  group_by(time) |>
  filter(row_number() == 1) |>
  ungroup() |>
  select(-location, -time, -nb_departure) %>% 
  as.matrix() ->
  out$X
#-----------------------------
out$fit_data <- list()
out$fit_data$W_valid <- utils$MC_train_test_split(!is.na(Y), validp, seed)
out$fit_data$train <- utils$to_incomplete(out$Y * out$fit_data$W_valid)
out$fit_data$valid <- out$Y[out$fit_data$W_valid==0]
out$fit_data$Y <- utils$to_incomplete(out$Y)


print(which(apply(out$fit_data$train, 1, function(x) all(x==0))))
print(which(apply(out$fit_data$train, 2, function(x) all(x==0))))
print((sum(out$fit_data$W_valid==0)/length(out$W)) + 
        (sum(out$fit_data$train!=0)/length(out$W)))
print(sum(!is.na(out$Y))/length(out$W))

return(out)

#--------------------------------------------------------
}


