library(data.table)
library(ggplot2)
library(reshape2)
library(corrplot)
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
source("./code_files/import_lib.R")
source("./BIXI/data-raw/bixi_data.R")
#----
model.dat <- load_bixi_dat(transpose = F)$model
#-----------------
Y = t(as.matrix(model.dat$depart))
Y[Y==0] = NA
X <- model.dat$X


impute_na_with_mean <- function(mat) {
 for (i in 1:ncol(mat)) {
  mat[is.na(mat[, i]), i] <- mean(mat[, i], na.rm = TRUE)
 }
 return(mat)
}

Y_imputed <- impute_na_with_mean(Y)
Y_long <- as.data.table(melt(Y_imputed))
colnames(Y_long) <- c("Location", "Day", "Departures")
Y_long <- Y_long[, .(Location, Day = as.integer(Day), Departures)]

X_long <- as.data.table(melt(X))
colnames(X_long) <- c("Day", "Covariate", "Value")
X_long <- X_long[, .(Day = as.integer(Day), Covariate, Value)]

merged_data <- merge(Y_long, X_long, by = "Day", allow.cartesian=TRUE)
correlations <- merged_data[, .(Correlation = cor(Departures, Value)), by = Covariate]
print(arrange(correlations,desc(abs(Correlation))))

cor_matrix <- dcast(merged_data, Location + Day ~ Covariate, value.var = "Value")
cor_matrix <- cor_matrix[, -c(1, 2)]
cor_matrix <- cor(cor_matrix, use = "complete.obs")

corrplot(cor_matrix, method = "circle")

cor_matrix[upper.tri(cor_matrix,T)] <- NA
round(cor_matrix,2) |> kable()

# ggplot(merged_data, aes(x = Value, y = Departures)) +
#  geom_point(alpha = 0.5) +
#  facet_wrap(~ Covariate, scales = "free") +
#  labs(title = "Scatter Plots of Covariates vs. Departures",
#       x = "Covariate Value",
#       y = "Number of Departures") +
#  theme_minimal()


#summary(Y_imputed)

summary(X)
X
#################
# Reshape Y and X for analysis
Y_long <- as.data.table(melt(Y_imputed))
colnames(Y_long) <- c("Location", "Day", "Departures")
Y_long <- Y_long[, .(Location, Day = as.integer(Day), Departures)]

X_long <- as.data.table(melt(X))
colnames(X_long) <- c("Day", "Covariate", "Value")
X_long <- X_long[, .(Day = as.integer(Day), Covariate, Value)]
merged_data <- merge(Y_long, X_long, by = "Day", allow.cartesian=TRUE)
average_departures <- merged_data[, .(Average_Departures = median(Departures)), by = .(Covariate, Value)]

ggplot(average_departures, aes(x = Value, y = Average_Departures)) +
 #geom_line() +
 geom_point() +
 facet_wrap(~ Covariate, scales = "free_x") +
 labs(title = "Median Number of Departures vs. Covariates",
      x = "Covariate Value",
      y = "Median Number of Departures") +
 theme_minimal()


ggplot(average_departures, aes(x = Value^2, y = Average_Departures)) +
  #geom_line() +
  geom_point() +
  facet_wrap(~ Covariate, scales = "free_x") +
  labs(title = "Median Number of Departures vs. Covariates (squared)",
       x = "Covariate Value",
       y = "Median Number of Departures") +
  theme_minimal()

average_departures %>%
  as.data.frame() |> 
  group_by(Covariate) |> 
  mutate(Value = (Value - median(Value))/IQR(Value)) |> 
  ungroup() |>

ggplot(aes(x = Value, y = Average_Departures)) +
  #geom_line() +
  geom_point() +
  facet_wrap(~ Covariate, scales = "free_x") +
  labs(title = "Median Number of Departures vs. Covariates (Normalized)",
       x = "Covariate Value",
       y = "Median Number of Departures") +
  theme_minimal()


average_departures %>%
  as.data.frame() |> 
  group_by(Covariate) |> 
  mutate(Value = (Value - median(Value))/IQR(Value)) |> 
  ungroup() |> 
  
  ggplot(aes(x = Value^2, y = Average_Departures)) +
  #geom_line() +
  geom_point() +
  facet_wrap(~ Covariate, scales = "free_x") +
  labs(title = "Median Number of Departures vs. Covariates (Normalized,squared)",
       x = "Covariate Value",
       y = "Median Number of Departures") +
  theme_minimal()
#---------------------------------------------------------------
med_scale <- function(x) (x - median(x)) / IQR(x)
X |> 
  as.data.frame() |> 
  transmute(num_metro_stations_log = log1p(num_metro_stations),
         num_restaurants_scaled = med_scale(num_restaurants),
         num_restaurants_sq_scaled = med_scale(num_restaurants^2),
         num_other_commercial_scaled = med_scale(num_other_commercial),
         num_pop_log = log1p(num_pop),
         num_bus_stations_scaled = med_scale(num_bus_stations),
         num_bus_routes_scaled = med_scale(num_bus_routes),
         num_bus_routes_sq_scaled = med_scale(num_bus_routes^2),
         num_university_bin = num_university,
         area_park_log = log1p(area_park),
         len_cycle_path_log = log1p(len_cycle_path),
         len_major_road_scaled = med_scale(len_major_road),
         len_minor_road_scaled = med_scale(len_minor_road),
         capacity_scaled = med_scale(capacity),
         walkscore_sq_scaled = med_scale(walkscore^2),
         walkscore_scaled = med_scale(walkscore))  |> 
  as.matrix() |> 
  melt() |> 
  as.data.table() -> 
  X_long
colnames(X_long) <- c("Day", "Covariate", "Value")
X_long <- X_long[, .(Day = as.integer(Day), Covariate, Value)]
merged_data <- merge(Y_long, X_long, by = "Day", allow.cartesian=TRUE)
average_departures <- merged_data[, .(Average_Departures = median(Departures)), by = .(Covariate, Value)]

ggplot(average_departures, aes(x = Value, y = Average_Departures)) +
  #geom_line() +
  geom_point() +
  facet_wrap(~ Covariate, scales = "free_x") +
  labs(title = "Median Number of Departures vs. Covariates",
       x = "Covariate Value",
       y = "Median Number of Departures") +
  theme_minimal()



correlations <- merged_data[, .(Correlation = cor(Departures, Value)), by = Covariate]
print(arrange(correlations,desc(abs(Correlation))))





###############
# Calculate daily total departures
daily_departures <- merged_data[, .(Total_Departures = sum(Departures)), by = Day]

# Plot the trend over time
ggplot(daily_departures, aes(x = Day, y = Total_Departures)) +
 geom_line() +
 labs(title = "Total Departures Over Time",
      x = "Day",
      y = "Total Departures") +
 theme_minimal()
#------------------------------------==

#---------------------------------------------------------------------
average_departures <- merged_data[, .(Average_Departures = mean(Departures)), by = .(Covariate, Value)]
ggplot(average_departures, aes(x = Value^2, y = Average_Departures)) +
 geom_point(alpha = 0.5) +
 geom_smooth(method = "lm", se = TRUE, color = "blue") +
 facet_wrap(~ Covariate, scales = "free_x") +
 labs(title = "Average Number of Departures vs. Covariates",
      x = "Covariate Value",
      y = "Average Number of Departures") +
 theme_minimal()

