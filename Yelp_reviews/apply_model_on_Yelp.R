
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

users <- readRDS("./Yelp_reviews/data/subset_PA/sample/users.RDS")
reviews <- readRDS("./Yelp_reviews/data/subset_PA/sample/reviews.RDS")
business <- readRDS("./Yelp_reviews/data/subset_PA/sample/business.RDS")

split_p <- list(train=0.6, test=0.2, valid=0.2)

n_rows <- nrow(reviews_matrix)
n_train <- floor(n_rows * train_proportion)

# Randomly select rows for the training set
set.seed(123)
train_indices <- sample(n_rows, n_train)

train_matrix <- reviews_matrix[train_indices, , drop = FALSE]
test_matrix <- reviews_matrix[-train_indices, , drop = FALSE]