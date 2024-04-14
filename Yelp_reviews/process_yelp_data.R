library(jsonlite)

file_path <- "./Yelp_reviews/data/yelp_academic_dataset_business.json"
businesses <- lapply(readLines(file_path, warn = FALSE), fromJSON)
business_df <- do.call(rbind, businesses) %>%  as.data.frame()
head(business_df)

business_subset <- business_df %>% 
 filter(state == "PA")
unlist(business_subset$review_count) %>%  sum()
business_to_retain = unlist(business_subset$business_id)

saveRDS(business_subset, "./Yelp_reviews/data/subset_PA/business_PA.RDS")

file_path <- "./Yelp_reviews/data/yelp_academic_dataset_review.json"
reviews <- lapply(readLines(file_path, warn = FALSE, n=10), fromJSON)
reviews <- do.call(rbind, reviews) %>%  as.data.frame()
head(reviews)

process_review <- function(line){
 review = fromJSON(line)
 if(review$business_id %in% business_to_retain){
  return(review)
 } else return(NULL)
}

matches <- list()
con <- file(file_path, open = "r")
on.exit(close(con))

while(TRUE) {
 line <- readLines(con, n = 1, warn = FALSE)
 if(length(line) == 0) break
 result <- process_review(line)
 if(!is.null(result)) {
  matches[[length(matches) + 1]] <- result
 }
}
reviews <- do.call(rbind, matches) %>%  as.data.frame()
head(reviews)

saveRDS(reviews, "./Yelp_reviews/data/subset_PA/reviews_PA.RDS")


reviews %>% dim()
head(reviews)
#-------------------------------
# extract user information the same way

users_to_retain = unique(unlist(reviews$user_id))
print(length(users_to_retain))


file_path <- "./Yelp_reviews/data/yelp_academic_dataset_user.json"
users <- lapply(readLines(file_path, warn = FALSE, n=10), fromJSON)
users <- do.call(rbind, users) %>%  as.data.frame()
head(users)


process_user <- function(line){
  user = fromJSON(line)
  if(user$user_id %in% users_to_retain){
    return(user)
  } else return(NULL)
}

matches <- list()
con <- file(file_path, open = "r")
on.exit(close(con))

while(TRUE) {
  line <- readLines(con, n = 1, warn = FALSE)
  if(length(line) == 0) break
  result <- process_user(line)
  if(!is.null(result)) {
    matches[[length(matches) + 1]] <- result 
  }
}
users <- do.call(rbind, matches) %>%  as.data.frame()
head(users)

saveRDS(users, "./Yelp_reviews/data/subset_PA/users_PA.RDS")
users %>% dim()
head(users)
#---------------------------------------------------------------------
# done writing. now read them back
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

users <- readRDS("./Yelp_reviews/data/subset_PA/users_PA.RDS")
reviews <- readRDS("./Yelp_reviews/data/subset_PA/reviews_PA.RDS")
business <- readRDS("./Yelp_reviews/data/subset_PA/business_PA.RDS")
#---------------------------------------------------------------------

unlist_df <- function(df){
  for (col_name in names(df)) {
    num.rows = nrow(df)
    # Check if the column is a list
    if (is.list(df[[col_name]])) {
      unlisted = unlist(df[[col_name]])
      if(length(unlisted)==num.rows)
        df[[col_name]] <- unlisted
    }
  }
  df
}

business <- unlist_df(business)
reviews <- unlist_df(reviews)
users <- unlist_df(users)
#------------------------------------------------------------------------
reviews %>%
  group_by(user_id) %>% 
  filter(n() >= 5) %>% 
  ungroup() %>% 
  group_by(business_id) %>% 
  filter(n() >= 5) %>% 
  ungroup() ->
  reviews
#-------------------------------------------------------------------
users %>%
  filter(user_id %in% reviews$user_id) ->
  users
business %>% 
  filter(business_id %in% reviews$business_id) ->
  business
reviews %>%
  filter(business_id %in% business$business_id) %>% 
  filter(user_id %in% users$user_id) ->
  reviews

#------------------------------------------------------------------------
business$city %>% table() %>%  as.data.frame() %>% arrange(desc(Freq)) %>% head(10)
business %>% 
  filter(city == "Philadelphia") ->
  business_philly

reviews %>%
  filter(business_id %in% business_philly$business_id) ->
  reviews_philly

users %>%
  filter(user_id %in% reviews_philly$user_id) ->
  users_philly
#-------------------------------------------------------------------------


#-------------------------------------------------------------------
# part 1 sampling: 1000 most popular resturants
reviews_philly %>%
  group_by(business_id) %>% 
  summarise(num_reviews = n()) %>% 
  arrange(desc(num_reviews)) %>% 
  slice_head(n=1000) %>%
  dplyr::select(business_id) ->
  sampled_business
sampled_business <- sampled_business$business_id

# out of those business, we sample the most popular 1000 users
reviews_philly %>%
  filter(business_id %in% sampled_business) %>% 
  group_by(user_id) %>% 
  summarise(num_reviews = n()) %>% 
  arrange(desc(num_reviews)) %>% 
  slice_head(n=1000) %>%
  select(user_id) ->
  sampled_users
sampled_users <- sampled_users$user_id

# we now keep only those:

reviews_philly %>%
  filter(business_id %in% sampled_business,
         user_id %in% sampled_users) ->
  reviews_philly
users_philly %>%
  filter(user_id %in% sampled_users) ->
  users_philly
business_philly %>%
  filter(business_id %in% sampled_business) ->
  business_philly
#---------------------------------
# rename ids
user_id_map <- setNames(seq_along(users_philly$user_id), users_philly$user_id)
business_id_map <- setNames(seq_along(business_philly$business_id),business_philly$business_id)
map_ids <- function(ids, map) {
  unname(map[as.character(ids)])
}
reviews_philly$user_id <- map_ids(reviews_philly$user_id, user_id_map)
reviews_philly$business_id <- map_ids(reviews_philly$business_id, business_id_map)
business_philly$business_id <- map_ids(business_philly$business_id, business_id_map) 
users_philly$user_id <- map_ids(users_philly$user_id, user_id_map)
#------------------------------------------------------------

reviews_matrix <- sparseMatrix(
  i = reviews_philly$user_id,
  j = reviews_philly$business_id,
  x = reviews_philly$stars,
  dims = c(nrow(users_philly), nrow(business_philly))
)
reviews_matrix <- as(reviews_matrix, "Incomplete")
length(reviews_matrix@x)/(nrow(reviews_matrix)*ncol(reviews_matrix))
#--------------------------------------------------------------------------
saveRDS(reviews_matrix, "Yelp_reviews/data/subset_PA/sample/reviews.RDS")
saveRDS(users_philly, "Yelp_reviews/data/subset_PA/sample/users.RDS")
saveRDS(business_philly, "Yelp_reviews/data/subset_PA/sample/business.RDS")
#---------------------------------------------------------------------------------
# end ......

#----------------------------------------------------------
# split train and tes
train_proportion <- 0.2
n_rows <- nrow(reviews_matrix)
n_train <- floor(n_rows * train_proportion)

# Randomly select rows for the training set
set.seed(123)
train_indices <- sample(n_rows, n_train)

train_matrix <- reviews_matrix[train_indices, , drop = FALSE]
test_matrix <- reviews_matrix[-train_indices, , drop = FALSE]

length(test_matrix@x) / length(reviews_matrix@x)
#-------------------------------------------------------------------------

business_df <- unlist_df(business_df)
business_df$city %>% table %>% as.data.frame %>% arrange(desc(Freq)) %>% head(10)
#-----------------------------------------------------------
# take a sample 


#-----------------------------------------------------------------
# user information
user.info <- users_philly %>% 
  dplyr::select(average_stars, review_count, useful, funny, cool, fans, contains("compliment")) %>% 
  as.matrix()

X_r = reduced_hat_decomp(user.info, tol=0.02)




sum(reviews_matrix==0)
unique(reviews_matrix@x)

dim(reviews_matrix)
max(reviews_philly$user_id, na.rm = T)
users_philly %>%
  filter(is.na(user_id))

head




