library(BKTR)
source("./code_files/import_lib.R")
# source("./BIXI/data-raw/bixi_data.R")

#' Generate BIXI data splits
#'
#' @param time_cov Logical; if TRUE, include time-varying covariates.
#' @param seed     Integer; random seed.
#' @param note     Character; appended to output file names.
#' @return          Saves `train.rds` and `test.rds` in `./BIXI/data/splits/`.
generate_bixi_data <- function(time_cov = TRUE, seed = 0, note = "") {
  # Set seed for reproducibility ---------------------------
  set.seed(seed)
  
  # Load raw data ------------------------------------------
  bixi_data  <- BixiData$new()
  data_df    <- bixi_data$data_df
  
  # Remove rows/locations with all-missing departures --------
  data_df <- data_df |>
    group_by(time) |>
    filter(!all(is.na(nb_departure))) |>
    ungroup() |>
    group_by(location) |>
    filter(!all(is.na(nb_departure))) |>
    ungroup()
  
  # Select covariates & reshape for matrix input -----------
  if (time_cov) {
    data_df <- data_df |>
      select(
        location,
        time,
        nb_departure,
        mean_temp_c,
        total_precip_mm,
        holiday
      ) |>
      arrange(location, time) |>
      rename(
        rows    = time,
        columns = location
      )
  } else {
    data_df <- data_df |>
      select(
        location,
        time,
        nb_departure,
        walkscore,
        capacity,
        num_metro_stations,
        num_bus_routes,
        num_university,
        len_minor_road,
        num_pop,
        area_park
      ) |>
      arrange(location, time) |>
      rename(
        rows    = location,
        columns = time
      )
  }
  
#------------------------------------------------------------------
# create splits for each of them:
# set.seed(0)
# data.df <- bixi.dat$data_df
# total_obs <- sum(! is.na(data.df$nb_departure))
# obs_indic <- which(! is.na(data.df$nb_departure))
# test_indic <- obs_indic[sample(1:total_obs, round(total_obs*0.3),replace = FALSE)]
# print(length(test_indic))
#----------------------------------------------------------------------------------

  # Initialize train/test -----------------------------------
  train_df <- data_df
  test_df  <- data_df
  
  # Determine dimensions & thresholds -----------------------
  num_rows    <- length(unique(train_df$rows))
  num_columns <- length(unique(train_df$columns))
  min_obs     <- if (time_cov) 4 else 80
  
  # Identify columns with enough data ----------------------
  eligible_columns <- train_df |>
    group_by(columns) |>
    filter(sum(is.na(nb_departure)) < min_obs) |>
    ungroup() |>
    distinct(columns) |>
    pull(columns)
  
  # Sample half of them for masking -------------------------
  set.seed(seed)
  effective_columns <- sample(
    eligible_columns,
    size    = round(num_columns * 0.5),
    replace = FALSE
  )
  
  # Mask test_df outside effective_columns -----------------
  missing_rate <- 0.95
  test_df <- test_df |>
    mutate(
      nb_departure = if_else(
        columns %in% effective_columns,
        nb_departure,
        NA_real_
      )
    )
  
  # Within each chosen column, split train/test entries -----
  for (col in effective_columns) {
    missing_rows <- sample(
      unique(train_df$rows),
      size    = floor(missing_rate * num_rows),
      replace = FALSE
    )
    
    train_df <- train_df |>
      mutate(
        nb_departure = if_else(
          columns == col & rows %in% missing_rows,
          NA_real_,
          nb_departure
        )
      )
    
    test_df <- test_df |>
      mutate(
        nb_departure = if_else(
          columns == col & !rows %in% missing_rows,
          NA_real_,
          nb_departure
        )
      )
  }
  
  # Add extra 20% random missing in train_df ---------------
  non_na_idx <- which(
    !is.na(train_df$nb_departure) &
      !train_df$columns %in% effective_columns
  )
  num_to_na <- round(0.2 * length(non_na_idx))
  set.seed(seed)
  na_idx <- sample(non_na_idx, num_to_na)
  train_df$nb_departure[na_idx] <- NA_real_
  
  # Finalize test_df (drop NA departures) -----------------
  test_df <- test_df |>
    filter(!is.na(nb_departure))
  
  # Save splits --------------------------------------------
  file_prefix <- paste0(
    "./BIXI/data/splits/split_",
    if (time_cov) "T" else "L",
    note,
    "_"
  )
  saveRDS(train_df, file = paste0(file_prefix, "train.rds"))
  saveRDS(test_df,  file = paste0(file_prefix, "test.rds"))
}