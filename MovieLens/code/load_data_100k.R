load_movielens_100k <-
  function(subset = "b",
           remove_bad_movies = TRUE,
           validp = 0.2,
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
      # dtrain %<>%
      #   #mutate(movie_id = as.integer(factor(movie_id)))
      #   mutate(movie_id = ifelse(movie_id < bad_movies[1], movie_id, movie_id -
      #                              1)) %>%
      #   mutate(movie_id = ifelse(movie_id < bad_movies[2], movie_id, movie_id -
      #                              1))
      # 
      # dtest %<>%
      #   #mutate(movie_id = as.integer(factor(movie_id)))
      #   mutate(movie_id = ifelse(movie_id < bad_movies[1], movie_id, movie_id -
      #                              1)) %>%
      #   mutate(movie_id = ifelse(movie_id < bad_movies[2], movie_id, movie_id -
      #                              1))
      
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
    print(sum(is.na(out$Y)) / length(out$W)- 
            (sum(is.na(out$O)) / length(out$W)) )
    print(sum(is.na(out$O)) / length(out$W))
    print(sum(is.na(out$Y)) / length(out$W))
    
    #-----------------------------
    out$fit_data <- list()
    out$fit_data$W_valid <- utils$MC_train_test_split(!is.na(out$Y), validp, seed)
    out$fit_data$train <- utils$to_incomplete(out$Y * out$fit_data$W_valid)
    out$fit_data$valid <- out$Y[out$fit_data$W_valid==0]
    out$fit_data$Y <- utils$to_incomplete(out$Y)
    
    
    print(which(apply(out$fit_data$train, 1, function(x) all(x==0))))
    print(which(apply(out$fit_data$train, 2, function(x) all(x==0))))
    print((sum(out$fit_data$W_valid==0)/length(out$W)) + 
            (sum(out$fit_data$train!=0)/length(out$W)))
    print(sum(!is.na(out$Y))/length(out$W))
    
    #----------------------------------------
    return(out)
  }
