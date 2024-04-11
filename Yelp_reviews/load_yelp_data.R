library(jsonlite)

file_path <- "./Yelp_reviews/data/yelp_academic_dataset_business.json"
businesses <- lapply(readLines(file_path, warn = FALSE), fromJSON)
business_df <- do.call(rbind, businesses) %>%  as.data.frame()
head(business_df)

business_subset <- business_df %>% 
 filter(state == "PA")
unlist(business_subset$review_count) %>%  sum()
business_to_retain = unlist(business_subset$business_id)



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

saveRDS(reviews, "reviews.RDS")


reviews %>% dim()
head(reviews)
#-------------------------------
# extract user information the same way

users_to_retain = unlist(business_subset$business_id)



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

saveRDS(reviews, "reviews.RDS")


reviews %>% dim()
head(reviews)