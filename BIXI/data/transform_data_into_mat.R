setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
# library(BKTR)
# source("./code_files/import_lib.R")


bdat <- BixiData$new()
i = 2
model.dat <- list()

train.df <-
 readRDS(paste0("./BIXI/data/splits/split_", i, "_train.rds"))
test.df <-
 readRDS(paste0("./BIXI/data/splits/split_", i, "_test.rds"))

train.df |>
 as.data.frame() |>
 mutate(time = as.Date(time)) |>
 arrange(location, time) ->  train.df

X <- train.df |>
 group_by(time) |>
 filter(row_number() == 1) |>
 ungroup() |>
 select(-location, -time, -nb_departure)

model.dat$X <- as.matrix(X)

reshape2::dcast(train.df, time ~ location, value.var = "nb_departure") |>
 select(-time) |>
 as.matrix() ->
 Y

model.dat$Y <- Y

obs.mask <- (1*!is.na(Y)) %>%  
 {
  colnames(.) <- NULL
  .
 } %>%
 as.matrix() 
print(sum(obs.mask == 1) / length(obs.mask))

train.df |>
 select(time, location, nb_departure) |>
 merge(
  select(test.df, time, location, nb_departure),
  by = c("time", "location"),
  all.x = T
 ) |>
 as.data.frame() |>
 arrange(location, as.Date(time))  |>
 mutate(missing =  !(is.na(nb_departure.x) &
                      (!is.na(nb_departure.y)))) -> mixed

print(sum(1 - mixed$missing) / nrow(mixed))

reshape2::dcast(mixed, time ~ location, value.var = "missing") %>%
 select(-time) %>%
 {
  colnames(.) <- NULL
  .
 } %>%
 as.matrix() ->
 #as.integer() ->
 test.mask

test.mask = 1 * test.mask

print(sum(1 - test.mask) / length(test.mask))
obs.mask[1:5,1:5]
test.mask[1:5,1:5]

valid_mask <- matrix.split.train.test(obs.mask, testp = 0.2)


model.dat$masks <- list(tr_val = obs.mask, test = test.mask, valid = valid_mask)
#-----------------------
# get test matrix
train.df |>
 select(time, location, nb_departure) |>
 merge(
  select(test.df, time, location, nb_departure),
  by = c("time", "location"),
  all.x = T
 ) |>
 as.data.frame() |>
 arrange(location, as.Date(time))  |>
 select("time", "location", nb_departure.y) %>% 
reshape2::dcast(time ~ location, value.var = "nb_departure.y") |>
 select(-time) |>
 as.matrix() %>% 
 to_incomplete() ->
 test

model.dat$splits <- list(
 
train = to_incomplete(Y*valid_mask),
valid = to_incomplete(Y * (1-valid_mask)),         #Y[valid_mask==0],
test = test,
Y = to_incomplete(Y)
)
print(length(model.dat$splits$train@x)/length(Y))
print(length(model.dat$splits$test@x))
print(model.dat$splits$valid@x %>% length)
#-------------------------------
# get all observed matrix
train.df |>
 select(time, location, nb_departure) |>
 merge(
  select(test.df, time, location, nb_departure),
  by = c("time", "location"),
  all.x = T
 ) |>
 as.data.frame() |>
 arrange(location, as.Date(time))  |>
 select("time", "location", nb_departure.x, nb_departure.y) %>%
 mutate(nb_departure = ifelse(is.na(nb_departure.x), nb_departure.y, nb_departure.x)) %>% 
 reshape2::dcast(time ~ location, value.var = "nb_departure") |>
 select(-time) |>
 as.matrix() %>% 
 to_incomplete() ->
 observed

print(length(observed@x))
model.dat$depart <- observed


#--------------------------------------------------------------------------------------------
CASMC_0_Bixi_Wrapper(model.dat)
CASMC_2_Bixi_Wrapper(model.dat)
CASMC_3a_Bixi_Wrapper(model.dat)
