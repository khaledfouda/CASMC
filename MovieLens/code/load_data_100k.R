load_movielens_100k <-
  function(subset = "a",
           remove_bad_movies = TRUE,
           scale = TRUE,
           seed = 2023)
  {
    #subset ="a"; remove_bad_movies = F; scale=T; seed=2023
    out <- list()
    set.seed(seed)
    dtrain <- read_tsv(
      paste0("MovieLens/data/ml-100k/u", subset, ".base"),
      col_names = c("user_id", "movie_id", "rating", "timestamp"),
      show_col_types = FALSE
    ) %>% 
      rename(row=user_id, column=movie_id) %>% 
      arrange(row, column)
    
    dtest <- read_tsv(
      paste0("MovieLens/data/ml-100k/u", subset, ".test"),
      col_names = c("user_id", "movie_id", "rating", "timestamp"),
      show_col_types = FALSE
    ) %>% 
      rename(row=user_id, column=movie_id) %>% 
      arrange(row, column)
    
    out$X <-
      read_delim(
      "MovieLens/data/ml-100k/u.user",
      delim = "|",
      col_names = c("user_id", "age", "sex", "occupation", "zipcode"),
      show_col_types = FALSE
    ) %>% 
      rename(row=user_id) %>% 
      arrange(row) %>% 
      select(-row, -zipcode) %>%
      mutate(age.sq = age ** 2) %>%
      fastDummies::dummy_cols(
        remove_first_dummy  = TRUE,
        ignore_na = T,
        remove_selected_columns = T
      ) %>%
      as.matrix() %>% 
      utils$scalers$minmax()   

    head(dtrain)
    dim(dtrain)
    dim(dtest)
    # transmute(age.sq = age ** 2,
    #           Male = sex == "M") %>%
    # # mutate(Male = (Male - mean(Male)) / sd(Male),
    # #        age = (age - mean(age)) / sd(age),
    # #        age.sq = (age.sq - mean(age.sq)) / sd(age.sq)) %>%
    #  select(-age.sq) %>%
    # as.matrix
    #
    
    #-------------------------- drop movies with no ratings
    if (remove_bad_movies) {
      bad_movies <- (dtest %>%
                       filter(!column %in% dtrain$column))$column
      print(bad_movies)
      dtest %<>% filter(!column %in% bad_movies)
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
    #---------------------------------------------------------------
    combined <- 
      rbind(mutate(dtrain, train=1),
            mutate(dtest, train=0)) %>% 
      mutate(rating = utils$scalers$standard(rating))
    head(combined)
    summary(combined$rating)
    sum(combined$rating==0)
    
    out$O <-
    combined %>% 
      arrange(row, column) %>% 
      reshape2::dcast(row~column, value.var="rating") %>% 
    select(-row) %>% 
      as.matrix() 
    
    out$Y <-
      combined %>%
      arrange(row, column) %>% 
      mutate(rating=ifelse(train==1, rating, NA)) %>% 
      reshape2::dcast(row~column, value.var="rating") %>% 
      select(-row) %>% 
      as.matrix() 
    
    out$Y[1:5,1:5]
    
    
    out$W <- !(is.na(out$Y) & (!is.na(out$O)))
    # out$W <-
    #   combined %>%
    #   arrange(row, column) %>% 
    #   mutate(rating=ifelse(train==0, rating, NA)) %>% 
    #   reshape2::dcast(row~column, value.var="rating") %>% 
    #   select(-row) %>% 
    #   as.matrix() 
    # out$W <- ifelse(is.na(out$W), 1, 0)
    out$W[1:5,1:5]
    
    print(sum(out$W==0) / length(out$W))
    print(sum(is.na(Y)) / length(out$W)-(sum(is.na(A)) / length(out$W)) )
    print(sum(is.na(A)) / length(out$W))
    print(sum(is.na(Y)) / length(out$W))
    
    #----------------------------------------
    # prepare sparse matrix and X decomposition
    user.id = as.numeric(factor(as.character(dtrain$user_id),
                                levels = sort(unique(dtrain$user_id))))
    movie.id = as.numeric(factor(as.character(dtrain$movie_id),
                                 levels = sort(unique(dtrain$movie_id))))
    
    glob_mean = mean(c(dtrain$rating, dtest$rating))
    glob_sd = sd(c(dtrain$rating, dtest$rating))
    dtrain$rating = (dtrain$rating - glob_mean) / glob_sd
    dtest$rating = (dtest$rating - glob_mean) / glob_sd
    
    dtrain.mat <- sparseMatrix(
      i = user.id,
      j = movie.id,
      x = dtrain$rating,
      dims = c(max(dtrain$user_id), max(dtrain$movie_id))
    )
    
    dtrain.mat <- as(dtrain.mat, "Incomplete")
    #---------------------------------------------------------
    # tmp: split test. 90% test set.
    
    # Wsplit = 0 -> W=1 = 1.1%, Wsplit = 1 & w=0: 94%; Wsplit=1 & W=1 = 4%
    # W==1: 5%.  test size: 0.5%
    #------------------------------------------------------------
    X_r <- reduced_hat_decomp(dX, 1e-2)
    #--------------------------------------------
    if (scale) {
      biScale.out <- biScaleMatrix(as.matrix(dtrain.mat))
      dtrain.Y = biScale.out$scaledMat
      dtrain.mat <- as(dtrain.Y, "Incomplete")
    } else{
      biScale.out = NULL
      dtrain.Y = as.matrix(dtrain.mat)
    }
    
    dW <- !is.na(dtrain.Y)
    sum(dW == 0) / length(dW)
    dtrain.Y[is.na(dtrain.Y)] = 0
    W_valid <- matrix.split.train.test(dW, testp = 0.2)
    Y_train = (dtrain.Y * W_valid)
    Y_train_inc = Y_train
    Y_train_inc[Y_train_inc == 0] = NA
    Y_train_inc = as(Y_train_inc, "Incomplete")
    Y_valid = dtrain.Y[W_valid == 0]
    #---------------------------------------------
    # X test as sparse
    dtest.inc <- sparseMatrix(
      i = dtest$user_id,
      j = dtest$movie_id,
      x = dtest$rating,
      dims = c(max(dtest$user_id), max(dtest$movie_id))
    )
    
    dtest.inc <- as(dtest.inc, "Incomplete")
    #---------------------------------------------
    # prepare output
    out <- list()
    out$X_r = X_r
    out$W = dW
    out$glob_mean = glob_mean
    out$glob_sd = glob_sd
    out$biScale = biScale.out
    out$train.df = dtrain
    out$train.mat = dtrain.Y
    out$train.inc = dtrain.mat
    out$test.df = dtest
    out$test.inc = dtest.inc
    out$valid = list()
    out$valid$W = W_valid
    out$valid$train = Y_train
    out$valid$test = Y_valid
    out$valid$train.inc = Y_train_inc
    return(out)
  }
