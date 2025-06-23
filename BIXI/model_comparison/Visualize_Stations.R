source("./code_files/import_lib.R")
source("./BIXI/data/transform_data_into_mat_other_2.R")
source("./BIXI/model_comparison/fit_wrappers_bixi2.R")


dat <- load_model_bixi_dat3(time_cov = FALSE, seed = 2023,validp = .2) 


out0 <- SImpute_Sim_Wrapper(dat) 
hpar_r <- CAMC_Ridge_hparams 
hpar_r$beta$lambda.grid <- seq(10,7, length.out=20)
out1 <- CAMC_Ridge_Sim_Wrapper(dat, trace = F, hpar=hpar_r, return_fit = T, max_cores = 20)

hpar_n = CAMC_Nuclear_hparams
#hpar$beta$n.lambda = 80
out2 <- CAMC_Nuclear_Sim_Wrapper(dat, trace=F, hpar=hpar_n, return_fit = T)


hpar_l = CAMC_Lasso_hparams
hpar_l$beta$n.lambda = 60
hpar_l$beta$lambda.max = 1
out3 <- CAMC_Lasso_Sim_Wrapper(dat, trace=F, hpar = hpar_l, return_fit = T)

out4 <- Naive_Sim_Wrapper(dat)

out5 <- Mao_Sim_Wrapper(dat)
#-------------------------------------
train.df <- readRDS(paste0("./BIXI/data/splits/split_L_train.rds"))

bixi.dat <- BixiData$new()
bixi.dat$temporal_positions_df %<>% 
 filter(time %in% unique(train.df$columns))

#-----------------------------------
temporal_kernel = BKTR::KernelSE$new()
temporal_kernel$set_positions(bixi.dat$temporal_positions_df)
temporal_kernel$kernel_gen() %>%  as.matrix() -> temp_kern
temp_kern[1:10,1:10]
#----------------------------------
spatial_kernel = BKTR::KernelMatern$new(smoothness_factor = 3)
spatial_kernel$set_positions(bixi.dat$spatial_positions_df)
spatial_kernel$kernel_gen() %>%  as.matrix() -> spt_kern
spt_kern[1:10,1:10]
#---------------------------------------
hpar_n$laplacian$S.b <- hpar_l$laplacian$S.b <-
 hpar_r$laplacian$S.b <- temp_kern
hpar_n$laplacian$S.a <- hpar_l$laplacian$S.a <- spt_kern

hpar_n$laplacian$lambda.b = .58
hpar_n$laplacian$lambda.a = .8222
out9 <- CAMC_Nuclear_Sim_Wrapper(dat, trace=F, hpar=hpar_n, return_fit = T)
out9$results$model <- "CAMC-Nuclear+Temporal+Spatial"
#-----------------------------------------------
dat <- load_model_bixi_dat3(2023,.2) 
#------------------------------------------------
hpar_n$laplacian$S.a <- temp_kern
hpar_n$laplacian$S.b <- spt_kern
hpar_n$laplacian$lambda.a = .4556
hpar_n$laplacian$lambda.b = .8222
out9 <- CAMC_Nuclear_Sim_Wrapper(dat, trace=F, hpar=hpar_n, return_fit = T)
out9$results$model <- "CAMC-Nuclear+Temporal+Spatial"
#-------------------------
hpar_l$laplacian$S.a <- temp_kern
hpar_l$laplacian$S.b <- spt_kern
hpar_l$laplacian$lambda.a = .4111
hpar_l$laplacian$lambda.b = .7778
out10 <- CAMC_Lasso_Sim_Wrapper(dat, trace=F, hpar = hpar_l, return_fit = T)
out10$results$model <- "CAMC-Lasso+Temporal+Spatial"
#-------------------------
dat$X %>% dim()
dim(out10$fit$fit$beta)

Xbeta <- dat$X %*% out10$fit$fit$beta

summary
plot(1:nrow(dat$Y), apply(Xbeta, 1, mean))

bixi.dat$spatial_positions_df %>% arrange(location) %>%  head()

library(ggmap)



map.dat <- arrange(bixi.dat$spatial_positions_df,location) %>% 
                mutate(effect = apply(Xbeta, 2, mean)) # %>% 
        mutate(effect = out10$fit$fit$beta[3,])


map.dat %<>% mutate(effect = effect > quantile(map.dat$effect,.75))

# 
# montreal_map <- get_stamenmap(
#         bbox = c(left = -73.72, bottom = 45.42, right = -73.45, top = 45.62), 
#         zoom = 12, 
#         maptype = "terrain"
# )

#register_stadiamaps("c62baf14-c12d-4bac-ac8e-f027e861d229")
montreal_map <- get_stadiamap(
        bbox = c(left = -73.72, bottom = 45.42, right = -73.45, top = 45.62),
        zoom = 12, maptype = "stamen_toner_lite"
)

#montreal_map <- get_map(location = c(lon = -73.5673, lat = 45.5017), zoom = 12)
ggmap(montreal_map) +
        geom_point(data = map.dat,
                   aes(x = longitude, y = latitude, color = effect),
                   size = 1, alpha = 0.8) +
         # scale_color_gradient(low = "yellow", high = "red",
         #                      name = "Average Effect") +
        labs(title = "Locations in Montreal with Average Covariate Effect",
             x = "Longitude", y = "Latitude") +
        theme_minimal()


# 
# library(ggplot2)
# library(osmdata)
# library(ggspatial)
# 

# Define the bounding box for Montreal
bbox <- getbb("Montreal, Canada")

# Fetch the map using OpenStreetMap
montreal_map <- opq(bbox = bbox) %>%
        add_osm_feature(key = "highway") %>%
        osmdata_sf()

# Extract the bounding box coordinates for ggplot
bbox_coords <- as.numeric(st_bbox(montreal_map$osm_lines))

ggplot() +
        geom_sf(data = montreal_map$osm_lines, inherit.aes = FALSE, color = "gray") +
        geom_point(data = df, aes(x = longitude, y = latitude, color = average_effect), size = 3, alpha = 0.8) +
        scale_color_gradient(low = "yellow", high = "red", name = "Average Effect") +
        coord_sf(xlim = c(bbox_coords[1], bbox_coords[3]), ylim = c(bbox_coords[2], bbox_coords[4]), expand = FALSE) +
        labs(title = "Locations in Montreal with Average Covariate Effect",
             x = "Longitude", y = "Latitude") +
        theme_minimal()
