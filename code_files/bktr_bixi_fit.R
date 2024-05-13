library(BKTR)

bdat <- BixiData$new()
model.dat <- list()

bdat$data_df %>% head

bdat$spatial_positions_df %>% head

bdatm <- bdat$data_df

bdatm$nb_departure %>% summary()
locations <- unique(bdatm$location)
time <- unique(bdatm$time)

bdatm <- bdat$data_df %>%
 arrange(location, time) %>% 
 group_by(location) %>% 
 mutate(time_id = 1:length(time)) %>% 
 ungroup() %>% 
 group_by(time) %>% 
 mutate(location_id = 1:length(locations)) %>% 
 ungroup() %>% 
 #mutate(nb_departure = ifelse(is.na(nb_departure), 0, nb_departure))
 filter(!is.na(nb_departure))



depart.mat <- sparseMatrix(
 i = bdatm$location_id,
 j = bdatm$time_id,
 x = bdatm$nb_departure
)
depart.mat <- as.matrix(depart.mat)
summary(as.vector(depart.mat))
depart.mat[depart.mat==0] = NA
depart.mat <- as(depart.mat, "Incomplete")
model.dat$depart <- depart.mat
#---------------------------------------------------------
#location covariates
bdatm %>%
 group_by(location) %>% 
 mutate(across(everything(), ~ n_distinct(.) == 1 )) %>%
 ungroup() %>% 
 select(which(sapply(.,function(x) x[1]==TRUE))) %>% 
 names() ->
 location_covariates

bdatm %>%
 select(location_id,all_of(location_covariates)) %>% 
 group_by(location_id) %>% 
 filter(row_number()==1) %>% 
 ungroup() %>% 
 select(-location_id) -> location_covariates
model.dat$X_row <- location_covariates
#-----------------------------------------------------------
# time covariates

bdatm %>%
 group_by(time) %>% 
 mutate(across(everything(), ~ n_distinct(.) == 1 )) %>%
 ungroup() %>% 
 select(which(sapply(.,function(x) x[1]==TRUE))) %>% 
 names() ->
 time_covariates

bdatm %>%
 select(time_id,all_of(time_covariates)) %>% 
 group_by(time_id) %>% 
 filter(row_number()==1) %>% 
 ungroup() %>% 
 select(-time_id) -> time_covariates
model.dat$X_col <- time_covariates
#----------------------------------------------------------------
# constructing the similarty matrices
# Create a matrix of day differences
dates <- as.Date(unique(bdatm$time))
day_diff_matrix <- outer(dates, dates, 
                         Vectorize(function(x, y) abs(as.numeric(difftime(x, y, units = "days")))))
alpha <- 0.2
similarity_matrix <- exp(-alpha * day_diff_matrix)
similarity_matrix <- round(similarity_matrix, 5)

print(similarity_matrix)
sum(similarity_matrix==0)/length(similarity_matrix)
time_similarity <- similarity_matrix
model.dat$sim_col <- time_similarity
#-----------------------------------------------------------------------
# similarity matrix for the location
# Function to calculate Haversine distance between two points
haversine <- function(lon1, lat1, lon2, lat2) {
 R <- 6371  # Earth radius in kilometers
 lon1 <- lon1 * pi / 180
 lon2 <- lon2 * pi / 180
 lat1 <- lat1 * pi / 180
 lat2 <- lat2 * pi / 180
 
 delta_lon <- lon2 - lon1
 delta_lat <- lat2 - lat1
 
 a <- sin(delta_lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta_lon/2)^2
 c <- 2 * atan2(sqrt(a), sqrt(1 - a))
 d <- R * c  # Distance in kilometers
 
 return(d)
}
locations <- bdat$spatial_positions_df
distance_matrix <- matrix(0, nrow = nrow(locations), ncol = nrow(locations))

# Calculate distances for each pair of locations
for (i in 1:nrow(locations)) {
 for (j in 1:nrow(locations)) {
  distance_matrix[i, j] <- haversine(locations$longitude[i], locations$latitude[i],
                                     locations$longitude[j], locations$latitude[j])
 }
}

# Convert distance to similarity: here using an exponential decay
alpha = 0.2
similarity_matrix <- exp(-alpha * distance_matrix) 
summary(as.vector(similarity_matrix))
sum(similarity_matrix <= .256537)/length(similarity_matrix)
print(similarity_matrix)
dim(similarity_matrix)
location_similarity <- similarity_matrix
model.dat$sim_row <- location_similarity
#-------------------------------------------------------------------------
# Summarize to check if all values are constant within each group
#summarise(across(starts_with("constant"), all)) %>%
# Filter groups where all covariates are constant
#filter(across(everything(), identity))

mapply(summary,model.dat$X_col)
mapply(summary,model.dat$X_row)
#------------------------------------------------------------------------
# creating masks
masks <- list()
masks$obs <- as.matrix(model.dat$depart != 0)

split_p <- list(test = 0.4, valid = 0.2)


masks$test <- matrix.split.train.test(masks$obs, split_p$test)
masks$valid <-
 matrix.split.train.test(masks$obs * masks$test, split_p$valid)


(round(sum(masks$test == 0) / sum(masks$obs == 1), 3) == split_p$test)
(round(sum(masks$valid == 0) / sum(masks$obs == 1 &
                                    masks$test == 1), 3) == split_p$valid)
sum(masks$test)
#----------------------------------------------------------------------
# split the data
model.dat$splits <- list()
model.dat$splits$train = model.dat$depart * masks$test * masks$valid
model.dat$splits$test = model.dat$depart * (1 - masks$test)
model.dat$splits$valid = model.dat$depart * (1 - masks$valid)
length(model.dat$splits$train@x)
length(model.dat$splits$test@x)
length(model.dat$splits$valid@x)















