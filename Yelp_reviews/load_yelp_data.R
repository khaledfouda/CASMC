robust_scale <- function(x) {
  median_x <- median(x, na.rm = TRUE)
  iqr_x <- IQR(x, na.rm = TRUE)
  return((x - median_x) / iqr_x)
}

scale=FALSE; split_p = list(test = 0.2, valid = 0.2); covariates="rows"; seed=2023; subset="_4x3" 

load_Yelp_data <-
  function(scale = TRUE,
           split_p = list(test = 0.2, valid = 0.2),
           covariates = "rows",
           seed = NULL,
           subset = "_4x3") {
    if (!is.null(seed))
      set.seed(seed)
    stopifnot(covariates %in% c("rows", "columns"))
    reviews <-
      readRDS(paste0("./Yelp_reviews/data/subset_PA/sample/reviews",subset,".RDS")) #%>%  t()
    
    
    
    if(covariates == "rows"){
      
    
    users <- readRDS(paste0("./Yelp_reviews/data/subset_PA/sample/users",subset,".RDS")) %>% 
      as.data.frame() %>% 
      mutate(elite_count = sapply(strsplit(elite, ","), length)) %>% 
      dplyr::select(average_stars,
                    review_count,
                    useful,
                    funny,
                    cool,
                    fans,
                    elite_count,
                    contains("compliment"),
                    -compliment_profile ,
                    -compliment_cute,
                    -compliment_list) %>%
      mutate(across(everything(), ~ ifelse(is.na(.), mean(.,na.rm=T), .))) %>%
      mutate(across(everything(), robust_scale)) %>%
      #scale() %>%
      as.matrix()
    
    X <- users
    
    }else{
      
    business <-
      readRDS(paste0("./Yelp_reviews/data/subset_PA/sample/business",subset,".RDS")) 
    
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
    df1 <- select(df1,all_of(categories_to_keep$name))
    
    
    df2 = business %>% select(attributes)
    attributes_expanded <- df2 %>%
      mutate(row_id = row_number()) %>%
      unnest_wider(attributes) %>%
      select(-row_id)

    attributes_expanded %>% summarise_all(function(d) length(unique(d))) %>% 
      pivot_longer(cols=everything(), names_to="name", values_to = 'value') %>% 
      filter(value <= 5) %>% 
      select(name) -> attributes_to_keep
    
    attributes_expanded %>% 
      select(all_of(attributes_to_keep$name)) %>% 
      #mutate(across(everything(), factor)) %>% 
        fastDummies::dummy_cols(remove_first_dummy  = TRUE,ignore_na=T,
                              remove_selected_columns=T) %>%  
      select(where(~ mean(is.na(.)) <= 0.1)) %>% 
      mutate(across(everything(), ~ ifelse(is.na(.), mean(.,na.rm=T), .))) ->
      df2
    
    business %>%
      select(stars, review_count, is_open) %>% 
      cbind(df1) %>% 
      select(where(~ mean(is.na(.)) <= 0.1)) %>% 
      mutate(across(everything(), ~ ifelse(is.na(.)|is.nan(.), mean(.,na.rm=T), .))) %>%
      cbind(df2) %>% 
      select(where(~ length(table(.)) >1)) %>% 
      mutate(across(everything(), robust_scale)) %>%
      #scale() %>% 
      as.matrix() -> 
      business
    
    X <- business
    reviews = t(reviews)
    }
    #-----------------------------------------------------------------------------------
    
    
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
    #reviews@x <-  log1p(reviews@x)
    #reviews@x
    # normalize and remove outliers
    #normalize = function(x)
    #(x - min(x)) / (max(x)-min(x)) + .1 #sd(x)
    #reviews@x <- normalize(reviews@x)
    #-------------------------------------------------------------------------
    # scale
    # reviews@x <- ifelse(reviews@x>=4, 1, -1)
    
    #reviews@x[reviews@x > 7] = 7  
    glob_mean = mean(reviews@x)
    glob_sd = sd(reviews@x)
    reviews@x = robust_scale(reviews@x)
    # reviews@x = (reviews@x - glob_mean) / glob_sd
    if (scale) {
      biScale.out <- biScaleMatrix(as.matrix(reviews),min_max_scale = T,normalize = F)
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
    
    # users$elite_count <-
    # users %>%
    #   select(elite, elite_count) %>%  view
    # 
    # X_cov <- users %>%
    #   dplyr::select(average_stars,
    #                 review_count,
    #                 useful,
    #                 funny,
    #                 cool,
    #                 fans,
    #                 contains("compliment")) %>%
    #   scale() %>%
    #   as.matrix()
    
    # X_cov <- business %>% 
    #   dplyr::select(stars, review_count, is_open) %>% 
    #   scale() %>% 
    #   as.matrix()
    
    
    
    # require(stats)
    # pca_result <- prcomp(X_cov, scale. = TRUE)
    # summary(pca_result)
    # pca_scores <- pca_result$x[, 1:5]
    X_r <- reduced_hat_decomp(X)
    
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
      X_r = X_r,
      glob_mean = glob_mean,
      glob_sd = glob_sd
    )
    return(out)
  }
