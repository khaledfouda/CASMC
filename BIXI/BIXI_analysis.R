library(data.table)
library(ggplot2)
library(reshape2)
library(corrplot)
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
source("./code_files/import_lib.R")
source("./BIXI/data-raw/bixi_data.R")
#----
model.dat <- load_bixi_dat(transpose = T)$model
#-----------------
Y <- t(as.matrix(model.dat$depart))
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
print(correlations)

cor_matrix <- dcast(merged_data, Location + Day ~ Covariate, value.var = "Value")
cor_matrix <- cor_matrix[, -c(1, 2)]
cor_matrix <- cor(cor_matrix, use = "complete.obs")

corrplot(cor_matrix, method = "circle")


ggplot(merged_data, aes(x = Value, y = Departures)) +
 geom_point(alpha = 0.5) +
 facet_wrap(~ Covariate, scales = "free") +
 labs(title = "Scatter Plots of Covariates vs. Departures",
      x = "Covariate Value",
      y = "Number of Departures") +
 theme_minimal()


#summary(Y_imputed)

summary(X)


#################
# Reshape Y and X for analysis
Y_long <- as.data.table(melt(Y_imputed))
colnames(Y_long) <- c("Location", "Day", "Departures")
Y_long <- Y_long[, .(Location, Day = as.integer(Day), Departures)]

X_long <- as.data.table(melt(X))
colnames(X_long) <- c("Day", "Covariate", "Value")
X_long <- X_long[, .(Day = as.integer(Day), Covariate, Value)]

# Merge Y_long and X_long by Day
merged_data <- merge(Y_long, X_long, by = "Day", allow.cartesian=TRUE)

# Calculate average number of departures for each covariate level
average_departures <- merged_data[, .(Average_Departures = median(Departures)), by = .(Covariate, Value)]

ggplot(average_departures, aes(x = Value, y = Average_Departures)) +
 geom_line() +
 geom_point() +
 facet_wrap(~ Covariate, scales = "free_x") +
 labs(title = "Median Number of Departures vs. Covariates",
      x = "Covariate Value",
      y = "Median Number of Departures") +
 theme_minimal()


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
ggplot(average_departures, aes(x = Value, y = Average_Departures)) +
 geom_point(alpha = 0.5) +
 geom_smooth(method = "lm", se = TRUE, color = "blue") +
 facet_wrap(~ Covariate, scales = "free_x") +
 labs(title = "Average Number of Departures vs. Covariates",
      x = "Covariate Value",
      y = "Average Number of Departures") +
 theme_minimal()

