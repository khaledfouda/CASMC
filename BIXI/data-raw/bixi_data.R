library(data.table)


bixi_folder <- "./BIXI/"
bixi_spatial_features <- fread(paste0(bixi_folder,'data-raw/bixi_spatial_features.csv'), header = TRUE, encoding = 'UTF-8')
bixi_spatial_locations <- fread(paste0(bixi_folder,'data-raw/bixi_spatial_locations.csv'), header = TRUE, encoding = 'UTF-8')
bixi_station_departures <- fread(paste0(bixi_folder,'data-raw/bixi_station_departures.csv'), header = TRUE, encoding = 'UTF-8')
bixi_temporal_features <- fread(paste0(bixi_folder,'data-raw/bixi_temporal_features.csv'), header = TRUE, encoding = 'UTF-8')
bixi_temporal_locations <- fread(paste0(bixi_folder,'data-raw/bixi_temporal_locations.csv'), header = TRUE, encoding = 'UTF-8')

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
