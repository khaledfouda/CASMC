


setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
source("./code_files/import_lib.R")
source("./BIXI/data-raw/bixi_data.R")
#source("./BIXI/model_comparison/fit_wrappers_bixi2.R")

#generate_BIXI_data <- function(time_cov = TRUE, seed = 2023, note = ""){
        

bixi.dat <- BixiData$new()

# 
# bixi.dat$data_df[1:5,]
# 
# bktr.reg <- BKTRRegressor$new(
#  data_df = bixi.dat$data_df,
#  spatial_positions_df = bixi.dat$spatial_positions_df,
#  temporal_positions_df = bixi.dat$temporal_positions_df
# )
# 
# bktr.reg$mcmc_sampling()
# 
# 
# #saveRDS(bktr.reg, "./BIXI/results/BKTR_fit_backup.rds")
# 
# TSR$set_params(seed = 0, fp_type = "float32")
# 
# summary(bktr.reg)
# bktr.reg$beta_covariates_summary |> 
#  as.data.frame() |> 
#  arrange( desc(abs(Median)) )
# 

#---------------------------------
# divide into train/test:

set.seed(0)
data.df <- bixi.dat$data_df
total_obs <- sum(! is.na(data.df$nb_departure))
obs_indic <- which(! is.na(data.df$nb_departure))
test_indic <- obs_indic[sample(1:total_obs, round(total_obs*0.2),replace = FALSE)]
print(length(test_indic))
train.df <- data.df
train.df$nb_departure[test_indic] <- NA
test.df <- data.df[test_indic,]
# saveRDS(train.df, "split_train.rds")
# saveRDS(test.df, "split_test.rds")

data.df %<>%
        group_by(time) %>%
        filter(!all(is.na(nb_departure))) %>% 
        ungroup() %>% 
        group_by(location) %>%
        filter(!all(is.na(nb_departure))) %>% 
        ungroup()

#---------------------------------
data <- list()

data.df |> 
        as.data.frame() |> 
        select(location, time, nb_departure, walkscore, capacity, num_metro_stations,
               num_bus_routes, num_university, len_minor_road, num_pop, area_park) ->
        data[[1]]

data.df |> 
        as.data.frame() |> 
        select(location, time, nb_departure, mean_temp_c, total_precip_mm, holiday) ->
        data[[2]]

data.df |> 
        as.data.frame() |> 
        select(location, time, nb_departure, walkscore, capacity, num_metro_stations,
               num_bus_routes, num_university, len_minor_road, num_pop, area_park,
               mean_temp_c, total_precip_mm, holiday) ->
        data[[3]]
#------------------------------------------------------------------
# create splits for each of them:
# set.seed(0)
# data.df <- bixi.dat$data_df
# total_obs <- sum(! is.na(data.df$nb_departure))
# obs_indic <- which(! is.na(data.df$nb_departure))
# test_indic <- obs_indic[sample(1:total_obs, round(total_obs*0.3),replace = FALSE)]
# print(length(test_indic))
#----------------------------------------------------------------------------------
train <- test <- data[[2]]
head(train)
m = length(unique(train$time))
n = length(unique(train$location))


#effective_stations <- readRDS("./BIXI/data/splits/effective_stations.rds")
train %>%
        group_by(location) %>% 
        filter(sum(is.na(nb_departure)) <4 ) %>% 
        ungroup() %>% 
        select(location) %>% 
        unique() ->
        unique_locations

unique_locations$location -> unique_locations
length(unique_locations)

effective_stations  <- sample(unique_locations, round(n*.5), FALSE)
missing_rate = .95
test %<>% 
        mutate(nb_departure = replace(nb_departure,
                                      !location %in% effective_stations, NA))
for(station in effective_stations){
        missing_days = sample(unique(train$time), floor(missing_rate * m), FALSE)
        train %<>% 
                mutate(nb_departure = replace(nb_departure,
                                              location == station & time %in% missing_days, NA))
        test %<>% 
                mutate(nb_departure = replace(nb_departure,
                                              location == station & (!time %in% missing_days), NA))
        
}
sum(is.na(train$nb_departure)) / nrow(train)
sum(!is.na(test$nb_departure)) / nrow(train)
#---------------------------------------
train %>%
        group_by(location) %>%
        filter(all(is.na(nb_departure))) %>% 
        ungroup()
train %>%
        group_by(time) %>%
        filter(all(is.na(nb_departure))) %>% 
        ungroup()
# verify that no empty rows or columns
#----------------------------------------------------------------
#-------- for the rest, make extra (random) 20% missing but without adding them 
# to the test set.

non_na_indices <- which( (!is.na(train$nb_departure)) & 
                                 (!train$location %in% effective_stations) )
num_to_na <- round(0.2 * length(non_na_indices))
train$nb_departure[sample(non_na_indices, num_to_na, F)] <- NA
sum(is.na(train$nb_departure)) / nrow(train)
#------------------------------------------------------------

test %<>% filter(! is.na(nb_departure))
nrow(test) / nrow(train)

saveRDS(train,file =  paste0("./BIXI/data/splits/split_",1,"_train.rds"))
saveRDS(test,file =  paste0("./BIXI/data/splits/split_",1,"_test.rds"))

}
# END - IGNORE THE REST --- go to Desktop.
#-----------------------------------------------------------------------
for(i in 1:3){
        
        train.df <- data[[i]]
        train.df$nb_departure[test_indic] <- NA
        test.df <- data[[i]][test_indic,]
        saveRDS(train.df,file =  paste0("./BIXI/data/splits/split_",i,"_train.rds"))
        saveRDS(test.df,file =  paste0("./BIXI/data/splits/split_",i,"_test.rds"))
        
}
#-----------------------------------------------------------------------
for(i in 1:1){
        
        train.df <- readRDS(paste0("./BIXI/data/splits/split_",i,"_train.rds"))
        test.df <- readRDS(paste0("./BIXI/data/splits/split_",i,"_test.rds"))
        TSR$set_params(seed = 0, fp_type = "float32")
        train.df <- setkey(as.data.table(train.df), location, time)
        
        split1.fit <- BKTRRegressor$new(
                data_df = train.df,
                spatial_positions_df = bixi.dat$spatial_positions_df,
                temporal_positions_df = bixi.dat$temporal_positions_df
        )
        split1.fit$mcmc_sampling()
        saveRDS(split1.fit, paste0("./BIXI/data/fits/split_",i,"_fit.rds"))
        
}
#-----------------------------------------------------------------------
split1.fit$beta_covariates_summary

split1.fit |> summary()

split1.fit$imputed_y_estimates |> 
        as.data.frame() |> 
        merge(test.df, by = c("location", "time")) |> 
        select(location, time, y_est, nb_departure) -> 
        test.estimates


error_metric$rmse(test.estimates$y_est, test.estimates$nb_departure)
error_metric$mae(test.estimates$y_est, test.estimates$nb_departure)
error_metric$spearman(test.estimates$y_est, test.estimates$nb_departure)

