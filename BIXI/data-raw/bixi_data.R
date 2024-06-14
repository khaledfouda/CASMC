library(data.table)

generate_second_order_diff_similarities <- function(n) {
    Z <- diag(-2, n, n)
    diag(Z[-1,]) <- diag(Z[, -1]) <- 1
    Z[1, 1] <- Z[n, n] <-  -1
    Z
}

load_bixi_dat <-
    function(scale_response = F,
             scale_covariates = F,
             transpose = F,
             testp = 0.2,
             validp = 0.2,
             seed = NULL) {
        if(! is.null(seed)) set.seed(seed)
        bixi_folder <- "./BIXI/"
        bixi_spatial_features <-
            fread(
                paste0(bixi_folder, 'data-raw/bixi_spatial_features.csv'),
                header = TRUE,
                encoding = 'UTF-8'
            )
        bixi_spatial_locations <-
            fread(
                paste0(bixi_folder, 'data-raw/bixi_spatial_locations.csv'),
                header = TRUE,
                encoding = 'UTF-8'
            )
        bixi_station_departures <-
            fread(
                paste0(bixi_folder, 'data-raw/bixi_station_departures.csv'),
                header = TRUE,
                encoding = 'UTF-8'
            )
        bixi_temporal_features <-
            fread(
                paste0(bixi_folder, 'data-raw/bixi_temporal_features.csv'),
                header = TRUE,
                encoding = 'UTF-8'
            )
        bixi_temporal_locations <-
            fread(
                paste0(bixi_folder, 'data-raw/bixi_temporal_locations.csv'),
                header = TRUE,
                encoding = 'UTF-8'
            )
        
        setkey(bixi_spatial_features, location)
        setkey(bixi_spatial_locations, location)
        setkey(bixi_station_departures, location)
        setkey(bixi_temporal_features, time)
        setkey(bixi_temporal_locations, time)
        
        bixi.dat <- list(
            spatial_features = bixi_spatial_features,
            spatial_locations = bixi_spatial_locations,
            departures = bixi_station_departures,
            temporal_features = bixi_temporal_features,
            temporal_locations = bixi_temporal_locations
        )
        #---
        model.dat <- list()
        
        #bixi.dat$spatial_locations |> head()
        #bixi.dat$spatial_features |> head()
        #bdat$spatial_positions_df %>% head
        
        long_departures = melt((bixi.dat$departures),
                               id.vars = "location",
                               variable.name = "date",
                               value.name = "departures"
        ) |>
            as.data.frame() |>
            filter(!is.na(departures)) |>
            mutate(location = as.character(location),
                   date = as.Date(date)) |>
            arrange(location, date)
        
        
        #long_departures |> head()
        
        # location_factor <- factor(long_departures$location)
        # date_factor <- factor(long_departures$date)
        # locations <- as.character(unique(long_departures$location))
        # dates <- as.character(unique(long_departures$date))
        
        # sparse_mat <- sparseMatrix(
        #     i = as.integer(location_factor),
        #     j = as.integer(date_factor),
        #     x = long_departures$departures,
        #     dimnames = list(levels(location_factor), levels(date_factor))
        # )
        
        depart.mat <- pivot_wider(long_departures, 
                               names_from = location,
                               values_from = departures) |> 
            mutate(date = as.Date(date)) |> 
            arrange(date) 
        dates <- depart.mat$date
        depart.mat <-
            depart.mat |> 
            select(-date) |> 
            as.matrix() |> 
            t()
        locations <- rownames(depart.mat)
        colnames(depart.mat) <- as.character(dates) 
        summary(as.vector(depart.mat))
        
        depart.mat <- as(depart.mat, "Incomplete")
        model.dat$depart <- depart.mat
        #depart.mat[1:5,1:5]
        # length(sparse_mat@x)
        # dim(sparse_mat)
        # nrow(long_departures)
        # summary(sparse_mat@x)
        # length(unique(date_factor))
        
        #depart.mat <- as.matrix(sparse_mat)
        #depart.mat[depart.mat == 0] = NA
        #---------------------------------------------------------
        #location covariates
        # bixi.dat$spatial_features |> as.data.frame() |>
        #     head()
        med_scale <- function(x) x
        bixi.dat$spatial_features |>
            as.data.frame() |>
            mutate(location = as.character(location)) |> 
            filter(location %in% locations) |>
            mutate(order_vec = match(location, locations)) |> 
            arrange(order_vec) |>
            dplyr::select(-location, - order_vec) |>
        transmute(
                  #walkscore_sq_scaled = med_scale(walkscore^2),
                  walkscore_scaled = med_scale(walkscore),
                  len_minor_road_scaled = med_scale(len_minor_road),
                  num_restaurants_scaled = med_scale(num_restaurants),
                  #num_restaurants_sq_scaled = med_scale(num_restaurants^2),
                  capacity_scaled = med_scale(capacity),
                  area_park_log = log1p(area_park),
                  len_major_road_scaled = med_scale(len_major_road),
                  num_other_commercial_scaled = med_scale(num_other_commercial),
                  num_bus_stations_scaled = med_scale(num_bus_stations),
                  num_pop_log = log1p(num_pop),
                  num_university_bin = num_university,
                  num_metro_stations_log = log1p(num_metro_stations),
                  num_bus_routes_scaled = med_scale(num_bus_routes),
                  len_cycle_path_log = log1p(len_cycle_path)
                  #num_bus_routes_sq_scaled = med_scale(num_bus_routes^2)
                  ) |> 
            as.matrix() ->
            model.dat$X
        
        bixi.dat$spatial_features |>
            as.data.frame() |>
            mutate(location = as.character(location)) |> 
            filter(location %in% locations) |>
            mutate(order_vec = match(location, locations)) |> 
            arrange(order_vec) |>
            dplyr::select(-location, - order_vec) |>
            transmute(
                walkscore = (walkscore),
                len_minor_road = (len_minor_road),
                num_restaurants = (num_restaurants),
                capacity = (capacity),
                area_park = (area_park),
                len_major_road = (len_major_road),
                num_other_commercial = (num_other_commercial),
                num_bus_stations = (num_bus_stations),
                num_pop = (num_pop),
                num_university = num_university,
                num_metro_stations = (num_metro_stations),
                num_bus_routes = (num_bus_routes),
                len_cycle_path = (len_cycle_path)
            ) |> 
            as.matrix() ->
            model.dat$spatial_simple
        
        #-----------------------------------------------------------
        # time covariates
        # bixi.dat$temporal_features |> as.data.frame() |>
        #     head()
        
        bixi.dat$temporal_features |>
            as.data.frame() |>
            mutate(time = as.Date(time)) |>
            filter(time %in% dates) |>
            mutate(order_vec = match(time, dates)) |> 
            arrange(order_vec) |>
            dplyr::select(-time, -order_vec) |>
            as.matrix() ->
            model.dat$Z
        #-------------------------------------------------
        # similarity matrices
        model.dat$similarity.A <- NULL
        model.dat$similarity.B <-
            generate_second_order_diff_similarities(length(dates))
        
        #-------------------------------------------------
        if (scale_covariates) {
            model.dat$Z %<>% scale()
            model.dat$X %<>% scale()
        }
        if (transpose) {
            tmp <- model.dat$Z
            model.dat$Z <- model.dat$X
            model.dat$X <- tmp
            model.dat$depart %<>% t()
            tmp <- model.dat$similarity.B
            model.dat$similarity.B <- model.dat$similarity.A
            model.dat$similarity.A <- tmp
        }
        if (scale_response) {
            scaled_response <-
                scalers(model.dat$depart@x, "standard")
            #    (model.dat$depart@x - min(model.dat$depart@x)+1e-10) / 
            #    (max(model.dat$depart@x) - min(model.dat$depart@x))
                # (model.dat$depart@x - mean(model.dat$depart@x)) / sd(model.dat$depart@x)
            stopifnot(sum(scaled_response == 0) == 0)
            model.dat$depart@x <- scaled_response
        }
        #------------------------------------------------------------------------
        # creating masks
        
        masks <- list()
        masks$obs <- as.matrix(model.dat$depart != 0)
        masks$test <-
            matrix.split.train.test(masks$obs, testp)
        masks$valid <-
            matrix.split.train.test(masks$obs * masks$test, validp)
        model.dat$masks <-  masks
        
        
        # (round(sum(masks$test == 0) / sum(masks$obs == 1), 3) == split_p$test)
        # (round(sum(masks$valid == 0) / sum(masks$obs == 1 &  masks$test == 1), 3) == split_p$valid)
        # sum(masks$test)
        
        #----------------------------------------------------------------------
        # split the data
        
        model.dat$splits <- list()
        model.dat$splits$train = (model.dat$depart * masks$test * masks$valid) |> 
            as("dgCMatrix")
        model.dat$splits$test = (model.dat$depart * (1 - masks$test)) |> 
            as("dgCMatrix")
        model.dat$splits$valid = (model.dat$depart * (1 - masks$valid)) |> 
            as("dgCMatrix")
        model.dat$depart <- (as.matrix(model.dat$depart) * masks$test)
        model.dat$depart[model.dat$depart == 0] <- NA
        model.dat$depart %<>% as("Incomplete")
        
        length(model.dat$splits$train@x)
        length(model.dat$splits$test@x)
        length(model.dat$splits$valid@x)
        length(model.dat$depart)
        model.dat$dates = dates
        model.dat$locations = locations
        return(list(model = model.dat, raw = bixi.dat))
    }
