source("./code_files/import_lib.R")
source("./BIXI/data/transform_data_into_mat_other_2.R")
source("./BIXI/model_comparison/fit_wrappers_bixi2.R")


dat <- load_model_bixi_dat3(time_cov = FALSE, seed = 2023,validp = .2) 


out0 <- SImpute_Sim_Wrapper(dat) 
hpar_r <- CASMC_Ridge_hparams 
hpar_r$beta$lambda.grid <- seq(10,7, length.out=20)
out1 <- CASMC_Ridge_Sim_Wrapper(dat, trace = F, hpar=hpar_r, return_fit = T, max_cores = 20)

hpar_n = CASMC_Nuclear_hparams
#hpar$beta$n.lambda = 80
out2 <- CASMC_Nuclear_Sim_Wrapper(dat, trace=F, hpar=hpar_n, return_fit = T)


hpar_l = CASMC_Lasso_hparams
hpar_l$beta$n.lambda = 60
hpar_l$beta$lambda.max = 1
out3 <- CASMC_Lasso_Sim_Wrapper(dat, trace=F, hpar = hpar_l, return_fit = T)

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
out9 <- CASMC_Nuclear_Sim_Wrapper(dat, trace=F, hpar=hpar_n, return_fit = T)
out9$results$model <- "CASMC-Nuclear+Temporal+Spatial"

#-------------------------
dat$X %>% dim()

Xbeta <- dat$X %*% out9$fit$fit$beta

plot(1:587, apply(Xbeta, 1, mean))
bixi.dat$spatial_positions_df %>% arrange(location) %>%  head()

library(sf)
library(ggmap)


map.dat <- arrange(bixi.dat$spatial_positions_df,location) %>% 
                 mutate(effect = apply(Xbeta, 1, mean))     
montreal_map <- get_map(location = c(lon = -73.5673, lat = 45.5017), zoom = 12)
ggmap(map.dat) +
        geom_point(data = df, aes(x = longitude, y = latitude, color = average_effect), size = 3, alpha = 0.8) +
        scale_color_gradient(low = "yellow", high = "red", name = "Average Effect") +
        labs(title = "Locations in Montreal with Average Covariate Effect",
             x = "Longitude", y = "Latitude") +
        theme_minimal()
