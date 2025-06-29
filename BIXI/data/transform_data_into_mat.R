
#setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
# source("./code_files/import_lib.R")

#' Load BIXI model data and build masks/splits
#'
#' @param time_cov Logical; if TRUE, use time-varying covariates
#' @return A list with matrices X, Y, masks, and splits.
load_model_bixi_dat <- function(time_cov = TRUE) {
  # Raw data object ------------------------------
  bdat <- BixiData$new()
  
  # Split identifier ------------------------------
  split_type <- if (time_cov) "T" else "S"
  
  # Read pre-saved train/test splits --------------
  train_df <- readRDS(paste0(
    "./BIXI/data/splits/split_",
    split_type,
    "_train.rds"
  ))
  test_df  <- readRDS(paste0(
    "./BIXI/data/splits/split_",
    split_type,
    "_test.rds"
  ))
  
  # Date conversion & sorting ---------------------
  if (time_cov) {
    train_df <- train_df %>%
      as.data.frame() %>% 
      rename(time = rows, location = columns) %>%
      mutate(time = as.Date(time)) %>%
      arrange(location, time)
    
    test_df <- test_df %>%
      as.data.frame() %>%
      rename(time = rows, location = columns) %>%
      mutate(time = as.Date(time)) %>%
      arrange(location, time)
    
  } else {
    train_df <- train_df %>%
      as.data.frame() %>%
      rename(time = columns, location = rows) %>%
      mutate(time = as.Date(time)) %>%
      rename(tmp_loc = time, time = location) %>%
      rename(location = tmp_loc) %>%
      arrange(location, time)
    
    test_df <- test_df %>%
      as.data.frame() %>%
      rename(time = columns, location = rows) %>%
      mutate(time = as.Date(time)) %>%
      rename(tmp_loc = time, time = location) %>%
      rename(location = tmp_loc) %>%
      arrange(location, time)
  }
  
  # Build covariate matrix X ----------------------
  X <- train_df %>%
    group_by(time) %>%
    slice(1) %>%
    ungroup() %>%
    select(-location, -time, -nb_departure)
  
  # Build response matrix Y -----------------------
  Y <- reshape2::dcast(
    train_df,
    time ~ location,
    value.var = "nb_departure"
  ) %>%
    select(-time) %>%
    as.matrix()
  
  # Observation mask -----------------------------
  obs_mask <- (!is.na(Y)) * 1
  print(sum(obs_mask == 1) / length(obs_mask))
  
  # Identify test entries via merge --------------
  mixed <- train_df %>%
    select(time, location, nb_departure) %>%
    merge(
      select(test_df, time, location, nb_departure),
      by = c("time", "location"),
      all.x = TRUE
    ) %>%
    as.data.frame() %>%
    arrange(location, time) %>%
    mutate(
      missing = !(is.na(nb_departure.x) & !is.na(nb_departure.y))
    )
  print(sum(1 - mixed$missing) / nrow(mixed))
  
  test_mask <- reshape2::dcast(
    mixed,
    time ~ location,
    value.var = "missing"
  ) %>%
    select(-time) %>%
    { colnames(.) <- NULL; . } %>%
    as.matrix() * 1
  print(sum(1 - test_mask) / length(test_mask))
  print(obs_mask[1:5, 1:5])
  print(test_mask[1:5, 1:5])
  
  # Validation mask via MC split -----------------
  valid_mask <- utils$MC_train_test_split(obs_mask, testp = 0.2)
  
  # Assemble core model.dat -----------------------
  model_dat <- list(
    X     = as.matrix(X),
    Y     = Y,
    masks = list(
      tr_val = obs_mask,
      test   = test_mask,
      valid  = valid_mask
    )
  )
  
  # Prepare test matrix ---------------------------
  test_mat <- train_df %>%
    select(time, location, nb_departure) %>%
    merge(
      select(test_df, time, location, nb_departure),
      by = c("time", "location"),
      all.x = TRUE
    ) %>%
    as.data.frame() %>%
    arrange(location, time) %>%
    select(time, location, nb_departure.y) %>%
    reshape2::dcast(
      time ~ location,
      value.var = "nb_departure.y"
    ) %>%
    select(-time) %>%
    as.matrix() %>%
    utils$to_incomplete()
  
  # Compile splits --------------------------------
  Xq <- qr(model_dat$X)
  model_dat$splits <- list(
    train = utils$to_incomplete(model_dat$Y * valid_mask),
    valid = utils$to_incomplete(model_dat$Y * (1 - valid_mask)),
    test  = test_mat,
    Y     = utils$to_incomplete(model_dat$Y),
    Xq      = qr.Q(Xq),
    Xr      = qr.R(Xq)
  )
  
  print(length(model_dat$splits$train@x) / length(model_dat$Y))
  print(length(model_dat$splits$test@x))
  print(length(model_dat$splits$valid@x))
  
  #-------------------------------
  # get all observed matrix
  # train_df %>%
  #   select(time, location, nb_departure) %>%
  #   merge(
  #     select(test_df, time, location, nb_departure),
  #     by = c("time", "location"),
  #     all.x = TRUE
  #   ) %>%
  #   as.data.frame() %>%
  #   arrange(location, time)  %>%
  #   select("time", "location", nb_departure.x, nb_departure.y) %>%
  #   mutate(nb_departure = ifelse(is.na(nb_departure.x), nb_departure.y, nb_departure.x)) %>%
  #   reshape2::dcast(time ~ location, value.var = "nb_departure") %>%
  #   select(-time) %>%
  #   as.matrix() %>%
  #   utils$to_incomplete() ->
  #   observed
  #
  # print(length(observed@x))
  # model_dat$depart <- Y#observed
  
  return(model_dat)
}

#--------------------------------------------------------------------------------------------
# usage examples:
# model_dat <- load_model_bixi_dat()
# SImpute_Bixi_Wrapper(model_dat)
# Mao_Bixi_Wrapper(model_dat)
# CAMC_0_Bixi_Wrapper(model_dat, train_on_all = TRUE)
# CAMC_2_Bixi_Wrapper(model_dat, train_on_all = TRUE)
# CAMC_3a_Bixi_Wrapper(model_dat, train_on_all = TRUE) -> results
# Naive_Bixi_Wrapper(model_dat)
