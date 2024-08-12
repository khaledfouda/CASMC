setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
source("./code_files/import_lib.R")
source("./BIXI/data-raw/bixi_data.R")
#source("./BIXI/model_comparison/fit_wrappers_bixi2.R")

generate_BIXI_data <- function(time_cov = TRUE, seed = 0, note = ""){
        
set.seed(seed)

bixi.dat <- BixiData$new()
#---------------------------------
# divide into train/test:

data.df <- bixi.dat$data_df
#total_obs <- sum(! is.na(data.df$nb_departure))
#obs_indic <- which(! is.na(data.df$nb_departure))
#test_indic <- obs_indic[sample(1:total_obs, round(total_obs*0.2),replace = FALSE)]
#print(length(test_indic))
#train.df <- data.df
#train.df$nb_departure[test_indic] <- NA
#test.df <- data.df[test_indic,]
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

if(time_cov){
        data.df %<>% 
                as.data.frame() %>%  
                select(location, time, nb_departure, mean_temp_c,
                       total_precip_mm, holiday) %>% 
                arrange(location, time) %>% 
                rename(rows = time, columns = location)
        
}else{
        data.df %<>% 
                as.data.frame() %>%  
                select(location, time, nb_departure, walkscore, 
                       capacity, num_metro_stations,
                       num_bus_routes, num_university, len_minor_road, 
                       num_pop, area_park) %>% 
                arrange(location, time) %>% 
                rename(rows = location, columns = time)

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

train <- test <- data.df
head(train)
m = length(unique(train$rows))
n = length(unique(train$columns))
min_obs <- ifelse(time_cov, 4,  80)

train %>%
        group_by(columns) %>% 
        #mutate(s = sum(is.na(nb_departure))) %>%  select(s) %>%  summary()
        filter(sum(is.na(nb_departure)) < min_obs ) %>% 
        ungroup() %>% 
        select(columns) %>% 
        unique() ->
        unique_columns

unique_columns$columns -> unique_columns
length(unique_columns)

effective_columns  <- sample(unique_columns, round(n*.5), FALSE)
missing_rate = .95
test %<>% 
        mutate(nb_departure = replace(nb_departure,
                                      !columns %in% effective_columns, NA))
for(col in effective_columns){
        missing_rows = sample(unique(train$rows), floor(missing_rate * m), FALSE)
        train %<>% 
                mutate(nb_departure = replace(nb_departure,
                                              columns == col & rows %in% missing_rows, NA))
        test %<>% 
                mutate(nb_departure = replace(nb_departure,
                                              columns == col & (!rows %in% missing_rows), NA))
        
}
train %>%
        group_by(rows) %>% 
        summarise(non_missing_in_rows = sum(! is.na(nb_departure))) %>% 
        select(non_missing_in_rows) %>% 
        summary() %>%  cbind(
train %>%
        group_by(columns) %>% 
        summarise(non_missing_in_cols = sum(! is.na(nb_departure))) %>% 
        select(non_missing_in_cols) %>% 
        summary()
        ) %>% print()

print(sum(is.na(train$nb_departure)) / nrow(train))
print(sum(!is.na(test$nb_departure)) / nrow(train))
#---------------------------------------
train %>%
        group_by(rows) %>%
        filter(all(is.na(nb_departure))) %>% 
        ungroup() %>% print()
train %>%
        group_by(columns) %>%
        filter(all(is.na(nb_departure))) %>% 
        ungroup() %>%  print()
# verify that no empty rows or columns
#----------------------------------------------------------------
#-------- for the rest, make extra (random) 20% missing but without adding them 
# to the test set.

non_na_indices <- which( (!is.na(train$nb_departure)) & 
                                 (!train$columns %in% effective_columns) )
num_to_na <- round(0.2 * length(non_na_indices))
train$nb_departure[sample(non_na_indices, num_to_na, F)] <- NA
sum(is.na(train$nb_departure)) / nrow(train)
#------------------------------------------------------------

test %<>% filter(! is.na(nb_departure))
print(nrow(test) / nrow(train))

file_name = paste0("./BIXI/data/splits/split_",ifelse(time_cov,"T","L"),note,
                   "_")
print(paste0("Files saved to ", file_name, "..."))
saveRDS(train,file =  paste0(file_name, "train.rds"))
saveRDS(test,file =   paste0(file_name,  "test.rds"))

}
# END - IGNORE THE REST --- go to Desktop.
#-----------------------------------------------------------------------
