# -----------------------------------------------------------------------------#
# Project:  Combining satellite-derived forest deadwood estimates and 
#           vegetation modelling to predict climate impacts on Norway 
#           spruce forests in Central Europe
#            
# Author:   Tim Anders
# Email:    tim.anders@senckenberg.de
# Date:     2026-01-07
# 
# Description: 
#           1. Model calibration
#           Calibration of the empirical GLM for
#           weather-related Norway spruce tree mortality
# -----------------------------------------------------------------------------#


#### 1. Load libraries #### 
library(terra)
library(ggplot2)
library(dplyr)
library(data.table)

#### 2. Define directories #### 
rm(list = ls())
wd <- "/storage/tanders/standing_deadwood_data/canheight_analysis/model_ger_3x3km/"
setwd(wd)

# -----------------------------------------------------------------------------#

#### 3. Load data ####
# Read data and combine in one final predictor variable set
predictors_2019 <- fread("./clim_pred_set_2019_3x3km_2102.csv")
predictors_2020 <- fread("./clim_pred_set_2020_3x3km_2102.csv")

predictors_2019$log_inc_deadwood_spruce <- log(predictors_2019$inc_deadwood_spruce)
predictors_2019$inc_deadwood_spruce <- round(predictors_2019$inc_deadwood_spruce, digits = 3)
predictors_2019 <- predictors_2019[inc_deadwood_spruce > 0]
predictors_2020$log_inc_deadwood_spruce <- log(predictors_2020$inc_deadwood_spruce)
predictors_2020$inc_deadwood_spruce <- round(predictors_2020$inc_deadwood_spruce, digits = 3)
predictors_2020 <- predictors_2020[inc_deadwood_spruce > 0]

# -----------------------------------------------------------------------------#

#### 4. Process training & testing data ####

# Add year column
predictors_2019$year <- 2019
predictors_2020$year <- 2020

# Reorder columns
setcolorder(predictors_2019, c("x", "y", "year", setdiff(names(predictors_2019), c("x", "y", "year"))))
setcolorder(predictors_2020, c("x", "y", "year", setdiff(names(predictors_2020), c("x", "y", "year"))))

# Remove temporary columns
predictors_2019[, c("log_inc_deadwood_spruce") := NULL]
predictors_2020[, c("log_inc_deadwood_spruce") := NULL]

predictors_2019_filt <- predictors_2019[spruce_cover_blick >= 0.1]
predictors_2020_filt <- predictors_2020[spruce_cover_blick >= 0.1]

predictors_2019_filt <-predictors_2019_filt[, c(1:4, 8:11, 15, 16, 26)]
predictors_2020_filt <-predictors_2020_filt[, c(1:4, 8:11, 15, 16, 26)]


# Function to sample 75% of the data (select randomly 75% of rows)
set.seed(8)
sample_train <- function(dt) {
  dt[sample(.N, size = 0.75 * .N)]
}

# Sample 75% from each year
train_2019_dt <- sample_train(predictors_2019_filt)
train_2020_dt <- sample_train(predictors_2020_filt)

# Check the number of rows
nrow(na.omit(train_2019_dt))  # 16958
nrow(na.omit(train_2020_dt))  # 16922
nrow(na.omit(predictors_2019_filt[!train_2019_dt, on = .(x, y)]))  # 5656
nrow(na.omit(predictors_2020_filt[!train_2020_dt, on = .(x, y)]))  # 5641

# Combine data of 2019 and 2020 directly
train_data <- rbindlist(list(train_2019_dt, train_2020_dt))
test_data <- rbindlist(list(
  predictors_2019_filt[!train_2019_dt, on = .(x, y)], 
  predictors_2020_filt[!train_2020_dt, on = .(x, y)]
))

head(train_data)
head(test_data)
nrow(na.omit(train_data)) # 33880
nrow(na.omit(test_data)) # 11297

all_data <- rbindlist(list(train_data, test_data))

# Create a function for the deviance explained
deviance_explained <- function(model) {
  null_deviance <- model$null.deviance
  residual_deviance <- model$deviance
  
  deviance_explained <- 1 - (residual_deviance / null_deviance)
  return(deviance_explained)
}

# Calibration dataset
calib_dt <- train_data[, .SD, .SDcols = !c("x", "y", "year")]
head(calib_dt)
nrow(calib_dt)

# Set the response variable
response_var <- "inc_deadwood_spruce"

# -----------------------------------------------------------------------------#

#### 5. Model calibration ####
calib_dt$tasmin_DJF <- calib_dt$tasmin_DJF - 273.15
calib_dt <- na.omit(calib_dt)
head(calib_dt)
nrow(calib_dt)

calib_filtered_df <- calib_dt[spruce_cover_blick >= 0.1]
nrow(calib_filtered_df)


# Correlation matrix
#Select only numeric columns
numeric_vars <- calib_filtered_df[, .SD, .SDcols = names(calib_filtered_df)[sapply(calib_filtered_df, is.numeric)]]

# Compute correlation matrix
cor_matrix <- cor(numeric_vars, use = "complete.obs", method = "pearson")

# View matrix
print(round(cor_matrix, 2))


### 5.1 Model 1 calibration ####
model_sat <- glm(formula = inc_deadwood_spruce ~  cwb_JJA_roll_min_5yrs
                 + rsds_anom_JJA_roll_mean_2yrs
                 + rsds_anom_MAM
                 + tasmin_DJF 
                 + tasmax_anom_JJA_1yr
                 + late_frost_days,
                 data = calib_filtered_df, 
                 family = Gamma(link = "log"))
model_sat

# Deviance explained
dev_expl <- deviance_explained(model_sat) # 0.4242582
dev_expl
