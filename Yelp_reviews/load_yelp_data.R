
load_Yelp_data <-
  function(scale = TRUE,
           split_p = list(test = 0.2, valid = 0.2),
           seed = NULL) {
    if (!is.null(seed))
      set.seed(seed)
    users <- readRDS("./Yelp_reviews/data/subset_PA/sample/users.RDS")
    reviews <-
      readRDS("./Yelp_reviews/data/subset_PA/sample/reviews.RDS") %>%  t()
    business <-
      readRDS("./Yelp_reviews/data/subset_PA/sample/business.RDS")
    
    mask_observ <- as.matrix(reviews != 0)
    
    
    split_p <- list(test = 0.2, valid = 0.2)
    mask_test <- matrix.split.train.test(mask_observ, split_p$test)
    mask_valid <-
      matrix.split.train.test(mask_observ * mask_test, split_p$valid)
    
    # (round(sum(mask_test == 0) / sum(mask_observ == 1), 3) == split_p$test)
    # (round(sum(mask_valid == 0) / sum(mask_observ == 1 &
    #                                    mask_test == 1), 3) == split_p$valid)
    # #----------------------------------------------------------------------------------
    # reviews@x %>%  range()
    # hist(reviews@x, breaks = 40)
    # summary(reviews@x)
    # sd(reviews@x)
    #---------------------------------------------------------------
    # normalize and remove outliers
    reviews@x[reviews@x > 5] = 5
    # normalize = function(x)
    #  (x - mean(x)) / sd(x)
    #reviews@x <- normalize(reviews@x)
    #-------------------------------------------------------------------------
    # scale
    if (scale) {
      biScale.out <- biScaleMatrix(as.matrix(reviews))
      reviews <- as(biScale.out$scaledMat, "Incomplete")
    } else{
      biScale.out = NULL
    }
    #-------------------------------------------------------------------------
    train_set = reviews * mask_test * mask_valid
    test_set = reviews * (1 - mask_test)
    valid_set = reviews * (1 - mask_valid)
    length(train_set@x)
    length(test_set@x)
    length(valid_set@x)
    #--------------------------------------------------------------------------
    X_cov <- users %>%
      dplyr::select(average_stars,
                    review_count,
                    useful,
                    funny,
                    cool,
                    fans,
                    contains("compliment")) %>%
      scale() %>%
      as.matrix()
    
    X_cov <- business %>% 
      dplyr::select(stars, review_count, is_open) %>% 
      scale() %>% 
      as.matrix()
    
    
    
    #require(stats)
    #pca_result <- prcomp(X_cov, scale. = TRUE)
    #summary(pca_result)
    #pca_scores <- pca_result$x[, 1:5]
    X_r <- reduced_hat_decomp(X_cov)
    
    #-----------------------------------------------------------------------------
    out = list(
      train.inc = reviews * mask_test,
      test.inc = test_set,
      W_obs = mask_observ,
      W_test = mask_test,
      valid = list(W = mask_valid,
                   train = train_set,
                   valid = valid_set),
      biScale = biScale.out,
      X_r = X_r
    )
    return(out)
  }
