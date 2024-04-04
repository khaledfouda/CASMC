load_movielens_100k <-
 function(subset = "a",
          remove_bad_movies = FALSE)
 {
  dtrain <- read_tsv(
   paste0("MovieLens/ml-100k/u", subset, ".base"),
   col_names = c("user_id", "movie_id", "rating", "timestamp"),
   show_col_types = FALSE
  )
  dtest <- read_tsv(
   paste0("MovieLens/ml-100k/u", subset, ".test"),
   col_names = c("user_id", "movie_id", "rating", "timestamp"),
   show_col_types = FALSE
  )
  
  dX <- read_delim(
   "MovieLens/ml-100k/u.user",
   delim = "|",
   col_names = c("user_id", "age", "sex", "occupation", "zipcode"),
   show_col_types = FALSE
  ) %>%
   arrange(user_id) %>%
   transmute(Male = sex == "M",
             (age - mean(age)) / sd(age)) %>%
   mutate(Male = (Male - mean(Male)) / sd(Male)) %>%
   as.matrix
  #-------------------------- drop movies with no ratings
  if (remove_bad_movies) {
   bad_movies <- (dtest %>%
                   filter(!movie_id %in% dtrain$movie_id))$movie_id
   bad_movies
   dtest %<>% filter(!movie_id %in% bad_movies)
   dtrain %<>% filter(!movie_id %in% bad_movies)
   # reset the index
   dtrain %<>%
    #mutate(movie_id = as.integer(factor(movie_id)))
    mutate(movie_id = ifelse(movie_id < bad_movies[1], movie_id, movie_id -
                              1)) %>%
    mutate(movie_id = ifelse(movie_id < bad_movies[2], movie_id, movie_id -
                              1))
   dtest %<>%
    #mutate(movie_id = as.integer(factor(movie_id)))
    mutate(movie_id = ifelse(movie_id < bad_movies[1], movie_id, movie_id -
                              1)) %>%
    mutate(movie_id = ifelse(movie_id < bad_movies[2], movie_id, movie_id -
                              1))
   
  }
  
  #----------------------------------------
  # prepare sparse matrix and X decomposition
  dtrain.mat <- sparseMatrix(
   i = dtrain$user_id,
   j = dtrain$movie_id,
   x = dtrain$rating,
   dims = c(max(dtrain$user_id), max(dtrain$movie_id))
  )
  
  dtrain.mat <- as(dtrain.mat, "Incomplete")
  X_r <- reduced_hat_decomp(dX, 1e-2)
  #--------------------------------------------
  biScale.out <- biScaleMatrix(as.matrix(dtrain.mat))
  dtrain.Y = biScale.out$scaledMat
  dtrain.mat <- as(dtrain.Y, "Incomplete")
  
  dW <- !is.na(dtrain.Y)
  sum(dW == 0) / length(dW)
  dtrain.Y[is.na(dtrain.Y)] = 0
  W_valid <- matrix.split.train.test(dW, testp = 0.2)
  Y_train = (dtrain.Y * W_valid)
  Y_valid = dtrain.Y[W_valid == 0]
  test_error <- RMSE_error
 #---------------------------------------------
 # prepare output
 out <- list() 
 out$X_r = X_r
 out$train_unscaled = dtrain
 out$W = dW
 out$biScale = biScale.out
 out$train.mat = dtrain.Y
 out$test = dtest
 out$train.inc = dtrain.mat
 out$valid = list()
 out$valid$W = W_valid
 out$valid$train = Y_train
 out$valid$test = Y_valid
 return(out)
 }










