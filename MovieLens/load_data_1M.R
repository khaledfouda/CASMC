load_movielens_1M <-
  function(scale = FALSE)
  {
    scale = FALSE
    ratings <-
      read.table(
        "MovieLens/ml-1m/ratings.dat",
        sep = ":",
        colClasses = c(NA, "NULL")
      ) %>%
      as.data.frame()
    names(ratings) <-
      c("user.id", "movie.id", "rating", "timestamp")
    ratings %<>%
      arrange(user.id, movie.id)
    #--------------------------------------------------------------------------
    users <-
      read.table("MovieLens/ml-1m/users.dat",
                 sep = ":",
                 colClasses = c(NA, "NULL")) %>%
      as.data.frame()
    names(users) <-
      c("user.id", "sex", "age", "occupation", "zip.code")
    users %<>%
      arrange(user.id) %>%
      select(-user.id, -zip.code) %>%
      mutate(age.sq = age ** 2) %>%
      fastDummies::dummy_columns(
        ignore_na = T,
        remove_selected_columns = T,
        remove_most_frequent_dummy = T
      ) %>%
      scale() %>%
      as.matrix()
    
    X_r <- reduced_hat_decomp(users, 1e-2)
    #---------------------------------------------------------------------------------------
    # prepare sparse matrix and X decomposition
    ratings.inc <- sparseMatrix(
      i = ratings$user.id,
      j = ratings$movie.id,
      x = scale(ratings$rating)[, 1],
      dims = c(max(ratings$user.id), max(ratings$movie.id))
    )  %>%
      as("Incomplete")
    #-----------------------------------------------------------------
    # head(users)
    # dim(users)
    # head(ratings)
    # length(ratings.inc@x)
    # dim(ratings.inc)
    
    
    #------------------------------------------------------------------------------------
    
    #---------------------------------------------------------
    # tmp: split test. 90% test set.
    
    # train = 60% test = 20% validation = 20%
    
    
    # Wsplit = 0 -> W=1 = 1.1%, Wsplit = 1 & w=0: 94%; Wsplit=1 & W=1 = 4%
    # W==1: 5%.  test size: 0.5%
    #------------------------------------------------------------
    
    #--------------------------------------------
    if (scale) {
      # do not scale (again) as it'll introduce new zeros.
      biScale.out <-
        biScaleMatrix(as.matrix(ratings.inc), min_max_scale = F)
      ratings.inc = biScale.out$scaledMat %>% as("Incomplete")
    } else{
      biScale.out = NULL
      #ratings.mat = as.matrix(ratings.inc)
    }
    
    W_missing = as.matrix(ratings.inc != 0) # 1 if observed (1M)
    W_test_valid <-
      matrix.split.train.test(W_missing, testp = 0.4) #  0 for test/valid and 1 for train/missing
    W_valid <-
      matrix.split.train.test((1 - W_test_valid) * W_missing , testp = 0.5) # 0 for valid and 1 for train/test/missing
    W_test <-
      1 -  (1 - W_test_valid) * (W_valid) # 0 for test and 1 for train/valid/missing
    
    # verification:
    # sum(W_test_valid == 0) / sum(W_missing == 1)
    # sum(W_test == 0) / sum(W_missing == 1)
    # sum(W_valid == 0) / sum(W_missing == 1)
    # sum(W_test == 0 & W_missing == 0) / sum(W_missing == 1)
    
    #-----------------------------------------------------------
    # prepare validation and test vectors:
    y_valid_vec = ratings.inc[W_valid == 0]
    y_test_vec = ratings.inc[W_test == 0]
    
    
    #---------------------------------------------------------
    # prepare training sets
    y = y_test = y_train = ratings.inc
    y[W_test == 0] = 0
    y_test[W_test == 1] = 0
    y_train[W_test_valid == 0] = 0
    #-------------------------------------------------------
    #---------------------------------------------
    # prepare output
    out <- list()
    out$X_r = X_r
    out$biScale = biScale.out
    
    out$masks = list()
    out$masks$missing = W_missing
    out$masks$test = W_test
    out$masks$valid = W_valid
    
    out$y = y
    out$y_test = y_test
    out$y_test_vec = y_test_vec
    
    out$valid = list()
    out$valid$train = y_train
    out$valid$valid = y_valid_vec
    return(out)
  }
