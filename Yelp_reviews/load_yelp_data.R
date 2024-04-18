
load_Yelp_data <-
  function(scale = TRUE,
           split_p = list(test = 0.2, valid = 0.2),
           covariates = "rows",
           seed = NULL) {
    if (!is.null(seed))
      set.seed(seed)
    stopifnot(covarites %in% c("rows", "columns"))
    reviews <-
      readRDS("./Yelp_reviews/data/subset_PA/sample/reviews.RDS") #%>%  t()
    
    
    users <- readRDS("./Yelp_reviews/data/subset_PA/sample/users.RDS") %>% 
      mutate(elite_count = sapply(strsplit(elite, ","), length)) %>% 
      dplyr::select(average_stars,
                    review_count,
                    useful,
                    funny,
                    cool,
                    fans,
                    contains("compliment")) %>%
      scale() %>%
      as.matrix()
      
    business <-
      readRDS("./Yelp_reviews/data/subset_PA/sample/business.RDS") 
    
    categs = strsplit(unlist(business$categories), ",")
  
    
    all_categories <- unique(unlist(strsplit(unlist(business$categories), ",\\s*")))
    all_categories[1:5]
    
    # Function to create indicators
    indicator_function <- function(category, categories) {
      as.integer(grepl(category, categories))
    }
    df1 = business %>% select(categories)
    # Apply the function to each category to create indicator columns
    for (category in all_categories) {
      df1[paste0("category_", gsub("[^[:alnum:]]", "", category))] <- mapply(indicator_function, category, df1$categories)
    }
    
    df1 %>% select(-categories) %>% summarise_all(sum) %>% 
      pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>% 
      arrange(desc(value)) %>%  
      filter(value > 100) %>% 
      select(name) -> categories_to_keep
    
    
    business %>%
      select(stars, review_count, is_open) %>% 
      cbind(select(df1,all_of(categories_to_keep$name)))
    
    # df1 %>%  %>%  cbind(business)
    df1 = business %>% select(attributes)
    
    
    attributes_expanded <- df1 %>%
      mutate(row_id = row_number()) %>% # Add an ID to preserve original row association
      unnest_wider(attributes) %>%
      select(-row_id) # Remove the row_id if you don't need it anymore
    
    head(attributes_expanded)
    
    attributes_expanded %>% summarise_all(function(d) length(unique(d))) %>% 
      pivot_longer(cols=everything(), names_to="name", values_to = 'value') %>% 
      filter(value <= 5) %>% 
      select(name) -> attributes_to_keep
    
    attributes_expanded %>% 
      select(all_of(attributes_to_keep$name)) %>% 
      mutate(across(everything(), factor)) %>% 
      
      view
    
    # If not already installed: install.packages("fastDummies")
    library(fastDummies)
    attributes_expanded <- fastDummies::dummy_cols(attributes_expanded, remove_first_dummy = TRUE)
    
    
    # Now 'attributes_expanded' has one column for each list element from the 'attributes' column
    # If necessary, you can then join this back with the original data frame 'df' by 'row_id'
    df <- df %>%
      mutate(row_id = row_number()) %>%
      left_join(attributes_expanded, by = "row_id") %>%
      select(-row_id)
    
    
    unlist(categs[1])[2]
  
  
  unlist(business$categories)[1]    
    
    if(covariates == "columns"){
      X = business
      reviews = t(reviews)
    }else
      X = users
    
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
    
    #reviews[reviews==0] <- NA
    #reviews@x[reviews@x > 5] = 5
    #reviews@x <-  log1p(reviews@x)
    #reviews@x
    # normalize and remove outliers
    #normalize = function(x)
    #(x - min(x)) / (max(x)-min(x)) + .1 #sd(x)
    #reviews@x <- normalize(reviews@x)
    #-------------------------------------------------------------------------
    # scale
    if (scale) {
      biScale.out <- biScaleMatrix(as.matrix(reviews),T,F)
      reviews <- as(biScale.out$scaledMat, "Incomplete")
    } else{
      biScale.out = NULL
    }
    
    
    
    #-------------------------------------------------------------------------
    train_set = reviews * mask_test * mask_valid
    test_set = reviews * (1 - mask_test)
    valid_set = reviews * (1 - mask_valid)
    # length(train_set@x)
    # length(test_set@x)
    # length(valid_set@x)
    #--------------------------------------------------------------------------
    
    users$elite_count <-
    users %>%
      select(elite, elite_count) %>%  view
    
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
    
    # X_cov <- business %>% 
    #   dplyr::select(stars, review_count, is_open) %>% 
    #   scale() %>% 
    #   as.matrix()
    
    
    
    require(stats)
    pca_result <- prcomp(X_cov, scale. = TRUE)
    summary(pca_result)
    pca_scores <- pca_result$x[, 1:5]
    X_r <- reduced_hat_decomp(pca_scores)
    
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
