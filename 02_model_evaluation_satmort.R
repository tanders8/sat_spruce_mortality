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
#     Evaluation of the GLM for weather-related Norway spruce tree 
#     mortality (VIP, correlation matrix, partial dependency plots)
# -----------------------------------------------------------------------------#


#### 1.  Load libraries ####
library(car)
library(terra)
library(ggplot2)
library(ggcorrplot)
library(dplyr) 
library(ggpubr)



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

# Correlation matrix
# 1. Remove 'spruce_cover_blick' from the data
numeric_vars <- calib_filtered_df[, !c("spruce_cover_blick"), with = FALSE]

# 2. Ensure all remaining columns are numeric (should already be the case)
numeric_vars <- numeric_vars[, lapply(.SD, as.numeric)]

# 3. Compute correlation matrix
cor_matrix <- cor(numeric_vars, use = "complete.obs", method = "pearson")


# Plot
cor_plot <- ggcorrplot(cor_matrix, 
                       method = "circle", 
                       type = "upper", 
                       lab = TRUE, 
                       lab_size = 3, 
                       colors = c("red", "white", "blue"),
                       title = paste0("Correlation matrix of training data \n(n=", nrow(calib_filtered_df),")"),
                       ggtheme = ggplot2::theme_minimal())

ggsave(filename = "correlation_matrix_plot.png", plot = cor_plot,
       width = 16, height = 16, units = "cm", dpi = 300)


## Prediction ###
# -----------------------------------------------------------------------------#
### Make predictions and plot mean value
predictors_1998_2020_dt <- fread("/storage/tanders/standing_deadwood_data/canheight_analysis/model_ger_3x3km/final_predictors_2001_2020_1402.txt", header = "auto")
head(predictors_1998_2020_dt)
setnames(predictors_1998_2020_dt, "cwb_roll_mean_GS_2yrs", "cwb_anom_roll_mean_GS_2yrs")
setnames(predictors_1998_2020_dt, "gpp_anom_JJA_roll_mean_2yrs", "gpp_JJA_anom_roll_mean_2yrs")
predictors_1998_2020_dt$tasmin_DJF <- predictors_1998_2020_dt$tasmin_DJF - 273.15

# Predict
predictors_1998_2020_dt[, prediction := predict(model_sat, newdata = predictors_1998_2020_dt, type = "response")]

# Summarize the mean prediction per year
# Add standard deviation to predictors_annual_mean
predictors_annual_mean <- predictors_1998_2020_dt[, 
                                                  .(mean_prediction = mean(prediction, na.rm = TRUE),
                                                    sd_prediction = sd(prediction, na.rm = TRUE)),  # Calculate standard deviation
                                                  by = year]

# Add standard deviation to observational_annual_mean
observational_annual_mean <- all_data[, 
                                      .(mean_observation = mean(inc_deadwood_spruce, na.rm = TRUE),
                                        sd_observation = sd(inc_deadwood_spruce, na.rm = TRUE)),  # Calculate standard deviation
                                      by = year]

# Merge both data frames by year if you need them in a single data frame
predictors_annual_mean <- merge(predictors_annual_mean, observational_annual_mean, by = "year", all = TRUE)
sub_pred_annual_mean <- predictors_annual_mean[19:20,]

rmse <- sqrt(mean((sub_pred_annual_mean[!is.na(mean_observation), mean_prediction - mean_observation])^2))
nrmse <- rmse / mean(sub_pred_annual_mean$mean_observation, na.rm = TRUE)
r_sqrt <- cor(sub_pred_annual_mean$mean_prediction, sub_pred_annual_mean$mean_observation, use = "complete.obs")^2
mae <- mean(abs(predictors_annual_mean[!is.na(mean_observation), mean_prediction - mean_observation]))


# Create a formatted string of predictor names and their coefficients
predictor_text <- paste0(
  paste(names(model_sat$coefficients), round(model_sat$coefficients, 4), sep = " = ", collapse = "\n")
)

p <- ggplot(predictors_annual_mean, aes(x = year)) +
  geom_line(aes(y = mean_prediction * 100, color = "Predicted"), size = 1) +
  geom_point(aes(y = mean_prediction * 100, color = "Predicted")) +
  geom_point(aes(y = mean_observation * 100, color = "Observed"), na.rm = TRUE) +  # Ignore NAs
  
  # Add standard deviation as error bars
  geom_errorbar(aes(ymin = (mean_prediction - sd_prediction) * 100, 
                    ymax = (mean_prediction + sd_prediction) * 100), 
                width = 0.2, color = "black") +  # Adjust the width as needed
  
  geom_errorbar(aes(ymin = (mean_observation - sd_observation) * 100, 
                    ymax = (mean_observation + sd_observation) * 100), 
                width = 0.2, color = "red") +  # Adjust the width as needed
  
  ylim(0, 20) +
  theme_minimal() +
  labs(title = "Mean predicted increase of spruce mortality from 2001 to 2020",
       x = "Year",
       y = "Mean spruce tree mortality increase [%]",
       color = "Legend") +  # Legend title
  
  # Annotate model performance metrics
  annotate("text", x = 2001, y = 17.5, label = paste0('dev. expl. = ', round(dev_expl, 4)), size = 5, hjust = 0) +
  annotate("text", x = 2001, y = 16.8, label = paste0("NRMSE = ", round(nrmse, 4)), hjust = 0, size = 5) +
  annotate("text", x = 2001, y = 16.1, label = paste0("MAE = ", round(mae, 4)), hjust = 0, size = 5) +
  
  # Add predictor names and coefficients in a text box
  annotate("text", x = 2010, y = 16.1, label = predictor_text, hjust = 0, size = 4, fontface = "italic") +
  
  scale_color_manual(values = c("Predicted" = "black", "Observed" = "red"))
p

# ------------------------------------------------------------------------------------------------------------ #
# Variable importance plots based on R2 (30 permutations)
# 
# Testing data
# ------------------------------------------------------------------------------------------------------------ #


# Transform into °C
test_data$tasmin_DJF <- test_data$tasmin_DJF - 273.15
test_data <- na.omit(test_data)
test_data_filtred <- test_data[spruce_cover_blick >= 0.1]
nrow(test_data_filtred)


# Get predicted probabilities from the original model
pred_model_test <- stats::predict(model_sat, newdata = test_data_filtred, type = "response")
test_data_filtred$pred <- pred_model_test

# Calculate Pearson's R² for the original model
r_squared_original <- cor(test_data_filtred$inc_deadwood_spruce, test_data_filtred$pred)^2
print(paste("Pearson's R² of original model:", r_squared_original))

# Function to calculate Pearson's R² after permuting a predictor multiple times
calculate_permuted_r_squared <- function(data, model_formula, predictor, test, nsim = 30) {
  r_squared_values <- numeric(nsim)  # Vector to store R² results for each permutation
  # Print the name of the current predictor being permuted
  print(paste("Permuting predictor:", predictor))
  
  for (i in 1:nsim) {
    # Print current iteration
    print(paste("Permutation iteration:", i, "of", nsim))
    
    # Create a permuted version of the dataset
    data_permuted <- data %>%
      mutate(!!sym(predictor) := sample(!!sym(predictor), n(), replace = FALSE))  # Permute the predictor
    
    # Fit the model on the permuted data
    model_permuted <- glm(model_formula, family = Gamma(link = "log"), data = data_permuted)
    
    # Get predicted probabilities from the permuted model
    pred_probs_permuted <- predict(model_permuted, newdata = test, type = "response")
    
    # Calculate Pearson's R² for the permuted model
    r_squared_values[i] <- cor(test$inc_deadwood_spruce, pred_probs_permuted)^2  # Calculate R²
  }
  
  return(r_squared_values)  # Return vector of R² values
}

# Define your model formula
model_formula <- as.formula("inc_deadwood_spruce ~ cwb_JJA_roll_min_5yrs + 
                             rsds_anom_JJA_roll_mean_2yrs + 
                             rsds_anom_MAM + 
                             tasmin_DJF + 
                             tasmax_anom_JJA_1yr + 
                             late_frost_days")

# List of predictors to permute
predictors <- c("cwb_JJA_roll_min_5yrs", "rsds_anom_JJA_roll_mean_2yrs", 
                "rsds_anom_MAM", "tasmin_DJF", 
                "tasmax_anom_JJA_1yr", "late_frost_days")

# Calculate R² for each permuted predictor 30 times
r_squared_permutations <- sapply(predictors, function(pred) {
  calculate_permuted_r_squared(calib_filtered_df, model_formula, pred, test_data_filtred, nsim = 30)
})

# Create a data frame to store results for boxplot
results_vip <- data.frame(
  Predictor = rep(predictors, each = 30),
  R_Squared_Permutation = as.vector(r_squared_permutations)
)

# Calculate the R² reduction for each permutation
results_vip <- results_vip %>%
  mutate(R_Squared_Reduction = r_squared_original - R_Squared_Permutation) %>%
  mutate(R_Squared_rel_Reduction = ((r_squared_original - R_Squared_Permutation) / r_squared_original) * 100)

# Reorder the predictors by the median R² reduction
results_vip <- results_vip %>%
  mutate(Predictor = reorder(Predictor, R_Squared_rel_Reduction, FUN = median))

# Create a data.frame
r_squared_permutations_df <- as.data.frame(r_squared_permutations)

# Check for significance
p_value <- colMeans(r_squared_permutations_df <= r_squared_original) # All predictors are highly important

# Compute mean R² drop for each predictor
mean_r2_permuted <- colMeans(r_squared_permutations_df, na.rm = TRUE)

# Compute relative reduction in R² (scaled importance score)
rel_importance <- (r_squared_original - mean_r2_permuted) / r_squared_original

# Print relative importance of predictors
print(rel_importance)

df_summary <- results_vip %>%
  group_by(Predictor) %>%
  summarise(
    R2_Perm_min = min(R_Squared_Permutation, na.rm = TRUE),
    R2_Perm_max = max(R_Squared_Permutation, na.rm = TRUE),
    R2_Perm_mean = mean(R_Squared_Permutation, na.rm = TRUE),
    
    R2_Red_min = min(R_Squared_Reduction, na.rm = TRUE),
    R2_Red_max = max(R_Squared_Reduction, na.rm = TRUE),
    R2_Red_mean = mean(R_Squared_Reduction, na.rm = TRUE),
    
    R2_relRed_min = min(R_Squared_rel_Reduction, na.rm = TRUE),
    R2_relRed_max = max(R_Squared_rel_Reduction, na.rm = TRUE),
    R2_relRed_mean = mean(R_Squared_rel_Reduction, na.rm = TRUE),
    .groups = "drop"
  )

print(df_summary)

# Visualize the results as boxplots, sorted by R² reduction
vip_test_gg <- ggplot(results_vip, aes(x = Predictor, y = R_Squared_rel_Reduction)) +
  geom_boxplot(fill = "lightgrey", outlier.colour = "black") +
  theme_minimal() +
  labs(title = "Variable importance: R² Reduction from 30 predictor variable permutations",
       x = "Predictors", y = "Reduction of R² [%]") +
  coord_flip() +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
vip_test_gg

ggsave("/storage/tanders/standing_deadwood_data/figures_manuscript_2/vip_test.jpeg", plot = vip_test_gg, 
       width = 24, height = 12, dpi = 500, units = "cm")
ggsave("/storage/tanders/standing_deadwood_data/figures_manuscript_2/FIG2.tif", plot = vip_test_gg, 
       width = 24, height = 12, dpi = 500, units = "cm")



# ------------------------------------------------------------------------------------------------------------ #
# Variable importance plots based on RMSE (30 permutations)
# 
# Testing data
# ------------------------------------------------------------------------------------------------------------ #

# Calculate RMSE for the original model
test_data$pred <- stats::predict(model_sat, newdata = test_data, type = "response")
rmse_original <- sqrt(mean((test_data$inc_deadwood_spruce - test_data$pred)^2))
print(paste("RMSE of original model:", round(rmse_original, 2)))

# Function to calculate RMSE after permuting a predictor multiple times
calculate_permuted_rmse <- function(data, model_formula, predictor, test, nsim = 30) {
  rmse_values <- numeric(nsim)  # Vector to store RMSE results for each permutation
  # Print the name of the current predictor being permuted
  print(paste("Permuting predictor:", predictor))
  
  for (i in 1:nsim) {
    # Print current iteration
    print(paste("Permutation iteration:", i, "of", nsim))
    
    # Create a permuted version of the dataset
    data_permuted <- data %>%
      mutate(!!sym(predictor) := sample(!!sym(predictor), n(), replace = FALSE))  # Permute the predictor
    
    # Fit the model on the permuted data
    model_permuted <- glm(model_formula, family = Gamma(link = "log"), data = data_permuted)
    
    # Get predicted probabilities from the permuted model
    pred_probs_permuted <- predict(model_permuted, newdata = test, type = "response")
    
    # Calculate RMSE for the permuted model
    rmse_values[i] <- sqrt(mean((test$inc_deadwood_spruce - pred_probs_permuted)^2))  # Calculate RMSE
  }
  
  return(rmse_values)  # Return vector of RMSE values
}

# Define your model formula
model_formula <- as.formula("inc_deadwood_spruce ~ cwb_JJA_roll_min_5yrs + 
                             rsds_anom_JJA_roll_mean_2yrs + 
                             rsds_anom_MAM + 
                             tasmin_DJF + 
                             tasmax_anom_JJA_1yr + 
                             late_frost_days")

# List of predictors to permute
predictors <- c("cwb_JJA_roll_min_5yrs", "rsds_anom_JJA_roll_mean_2yrs", 
                "rsds_anom_MAM", "tasmin_DJF", 
                "tasmax_anom_JJA_1yr", "late_frost_days")

# Calculate RMSE for each permuted predictor 30 times
rmse_permutations <- sapply(predictors, function(pred) {
  calculate_permuted_rmse(calib_dt, model_formula, pred, test_data, nsim = 30)
})

# Create a data frame to store results for boxplot
results_vip <- data.frame(
  Predictor = rep(predictors, each = 30),
  RMSE_Permutation = as.vector(rmse_permutations)
)

# Calculate the RMSE reduction for each permutation
results_vip <- results_vip %>%
  mutate(RMSE_Increase = RMSE_Permutation - rmse_original) %>%
  mutate(RMSE_rel_Increase = ((RMSE_Permutation - rmse_original) / rmse_original) * 100)

# Reorder the predictors by the median RMSE reduction
results_vip <- results_vip %>%
  mutate(Predictor = reorder(Predictor, RMSE_rel_Increase, FUN = median))

# Create a data.frame
rmse_permutations_df <- as.data.frame(rmse_permutations)

# Check for significance
p_value <- colMeans(rmse_permutations_df >= rmse_original) # Higher values indicate importance

# Compute mean RMSE for each permuted predictor
mean_rmse_permuted <- colMeans(rmse_permutations_df, na.rm = TRUE)

# Compute relative reduction in RMSE (scaled importance score)
rel_importance_rmse <- (rmse_original - mean_rmse_permuted) / rmse_original

# Print relative importance of predictors based on RMSE
print(rel_importance_rmse)

# Visualize the results as boxplots, sorted by RMSE reduction
vip_test_gg <- ggplot(results_vip, aes(x = Predictor, y = RMSE_rel_Increase)) +
  geom_boxplot(fill = "lightgrey", outlier.colour = "black") +
  theme_minimal(base_size = 14) +  # sets base text size to 14
  labs(
    title = "Variable Importance: RMSE Increase from 30 Predictor Variable Permutations",
    x = "Predictors",
    y = "Increase of RMSE [%]"
  ) +
  coord_flip() +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

vip_test_gg


ggsave("/storage/tanders/standing_deadwood_data/figures_manuscript_2/vip_test_rmse.jpeg", plot = vip_test_gg, 
       width = 24, height = 12, dpi = 500, units = "cm")
ggsave("/storage/tanders/standing_deadwood_data/figures_manuscript_2/FIGS3.tif", plot = vip_test_gg, 
       width = 24, height = 12, dpi = 500, units = "cm")

# ------------------------------------------------------------------------------------------------------------ #
# Variable importance plots based on R2 (30 permutations)
# 
# Training data
# ------------------------------------------------------------------------------------------------------------ #
# Based on training data

# Calculate R² for each permuted predictor 30 times
r_squared_permutations <- sapply(predictors, function(pred) {
  calculate_permuted_r_squared(calib_dt, model_formula, pred, calib_dt, nsim = 30)
})

# Create a data frame to store results for boxplot
results <- data.frame(
  Predictor = rep(predictors, each = 30),
  R_Squared_Permutation = as.vector(r_squared_permutations)
)

# Calculate the R² reduction for each permutation
results <- results %>%
  mutate(R_Squared_Reduction = r_squared_original - R_Squared_Permutation) %>%
  mutate(R_Squared_rel_Reduction = ((r_squared_original - R_Squared_Permutation) / r_squared_original) * 100)

# Reorder the predictors by the median R² reduction
results <- results %>%
  mutate(Predictor = reorder(Predictor, R_Squared_rel_Reduction, FUN = median))

# Create a data.frame
r_squared_permutations_df <- as.data.frame(r_squared_permutations)

# Check for significance
p_value <- colMeans(r_squared_permutations_df <= r_squared_original) # All predictors are highly important

# Compute mean R² drop for each predictor
mean_r2_permuted <- colMeans(r_squared_permutations_df, na.rm = TRUE)

# Compute relative reduction in R² (scaled importance score)
rel_importance <- (r_squared_original - mean_r2_permuted) / r_squared_original

# Print relative importance of predictors
print(rel_importance)


# Visualize the results as boxplots, sorted by R² reduction
vip_train_gg <- ggplot(results, aes(x = Predictor, y = R_Squared_rel_Reduction)) +
  geom_boxplot(fill = "lightgrey", outlier.colour = "black") +
  theme_minimal() +
  labs(title = "Variable Importance: R² Reduction from 30 Predictor Variable Permutations (training data)",
       x = "Predictors", y = "Reduction of R² [%]") +
  coord_flip()


ggsave("./figures_filtered/vip_train.jpeg", plot = vip_test_gg, width = 32, height = 16, units = "cm")

# Calculate the AUC for the original model
pred_model_1a <- stats::predict(model_1a, newdata = test, type = "response")
test$pred_model_1a <- pred_model_1a

# Use pROC for AUC calculation
roc_original <- PRROC::roc.curve(scores.class0 = test$pred_model_1a, weights.class0=test$dead,
                                 curve=F)  # Calculate ROC
auc_original <- roc_original$auc  # Calculate AUC
print(paste("AUC of original model:", auc_original))

# Function to calculate AUC after permuting a predictor
calculate_permuted_auc <- function(data, model_formula, predictor, test) {
  # Create a permuted version of the dataset
  data_permuted <- data %>%
    mutate(!!sym(predictor) := sample(!!sym(predictor), n()))  # Permute the predictor
  
  # Fit the model on the permuted data
  model_permuted <- glm(model_formula, family = binomial(link = "logit"), data = data_permuted)
  
  # Get predicted probabilities from the permuted model
  pred_probs_permuted <- predict(model_permuted, newdata = test, type = "response")
  
  test$pred_probs_permuted <- pred_probs_permuted
  
  # Calculate AUC for the permuted model using pROC
  roc_permuted <- PRROC::roc.curve(scores.class0 = test$pred_probs_permuted, weights.class0=test$dead,
                                   curve=F)  # Calculate ROC
  auc_permuted <- roc_permuted$auc  # Calculate AUC
  
  return(auc_permuted)  # Return the AUC of the permuted model
}

# Define your model formula
model_formula <- as.formula("dead ~ CWB_summer_dev_1yr + CWB_summer_dev_2yr + 
                             rsds_long_mean_1981_2010_summer + 
                             temp_mean_long_mean_1981_2010_winter + 
                             rsds_spring_dev_1yr + rsds_winter_dev_2yr + 
                             prec_winter_dev_0yr + temp_mean_autumn_dev_1yr + 
                             rsds_summer_dev_1yr")

# List of predictors to permute
predictors <- c("CWB_summer_dev_1yr", "CWB_summer_dev_2yr", 
                "rsds_long_mean_1981_2010_summer", "temp_mean_long_mean_1981_2010_winter", 
                "rsds_spring_dev_1yr", "rsds_winter_dev_2yr", 
                "prec_winter_dev_0yr", "temp_mean_autumn_dev_1yr", 
                "rsds_summer_dev_1yr")

# Calculate AUC for each permuted predictor
auc_permutations <- sapply(predictors, function(pred) {
  # Calculate the AUC for the model with the current predictor permuted
  calculate_permuted_auc(train, model_formula, pred, test)
})

# Create a data frame to compare original AUC with permuted AUCs
results <- data.frame(
  Predictor = predictors,
  AUC_Permutation = auc_permutations
)

# Add the original model AUC
results <- rbind(data.frame(Predictor = "Original Model", AUC_Permutation = auc_original), results)

# Calculate the reduction in AUC for each predictor
results <- results %>%
  mutate(AUC_Reduction = auc_original - AUC_Permutation) %>%
  arrange(desc(AUC_Reduction))  # Sort by AUC reduction

# Print the results
print(results)

# Optional: Visualize the results
ggplot(results[-1, ], aes(x = reorder(Predictor, -AUC_Reduction), y = AUC_Reduction)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Variable Importance: AUC Reduction from Permutation",
       x = "Predictors", y = "AUC Reduction") +
  coord_flip()


#### 30 permutations

#####################################################################################################################
### Partial dependency plots
#####################################################################################################################
# model sat
# 1 cwb_JJA_roll_min_5yrs                      
cwb_JJA_roll_min_5yrs <- pdp::partial(
  object = model_sat, 
  pred.var = "cwb_JJA_roll_min_5yrs", 
  train = calib_dt
)
cwb_JJA_roll_min_5yrs$yhat <- exp(cwb_JJA_roll_min_5yrs$yhat) * 100
plot_cwb_JJA_roll_min_5yrs <- autoplot(cwb_JJA_roll_min_5yrs) +
  theme_bw() +
  ylim(c(0, 30)) +
  ggtitle("cwb_JJA_roll_min_5yrs [mm]") +
  ylab("Mortality rate [%]") +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# 2 tasmin_DJF  
par_tasmin_DJF <- pdp::partial(
  object = model_sat, 
  pred.var = "tasmin_DJF", 
  train = calib_dt
)
par_tasmin_DJF$yhat <- exp(par_tasmin_DJF$yhat) * 100
plot_tasmin_DJF <- autoplot(par_tasmin_DJF) +
  theme_bw() +
  ylim(c(0, 30)) +
  ggtitle("tasmin_DJF [°C]") +
  ylab("Mortality rate [%]")+
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# 5 rsds_anom_JJA_roll_mean_2yrs     
par_rsds_anom_JJA_roll_mean_2yrs <- pdp::partial(
  object = model_sat, 
  pred.var = "rsds_anom_JJA_roll_mean_2yrs", 
  train = calib_dt
)
par_rsds_anom_JJA_roll_mean_2yrs$yhat <- exp(par_rsds_anom_JJA_roll_mean_2yrs$yhat) * 100
plot_rsds_anom_JJA_roll_mean_2yrs <- autoplot(par_rsds_anom_JJA_roll_mean_2yrs) +
  theme_bw() +
  ylim(c(0, 30)) +
  ggtitle(expression(rsds_anom_JJA_roll_mean_2yrs ~ "[" * Wm^{-2} * "]")) +
  ylab("Mortality rate [%]")+
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# 6 late_frosts
par_late_frosts <- pdp::partial(
  object = model_sat, 
  pred.var = "late_frost_days", 
  train = calib_dt
)
par_late_frosts$yhat <- exp(par_late_frosts$yhat) * 100
plot_late_frosts <- autoplot(par_late_frosts) +
  theme_bw() +
  ylim(c(0, 30)) +
  ggtitle("late frost days [days]") +
  ylab("Mortality rate [%]")+
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# 7 tasmax_anom_JJA_1yr       
par_tasmax_anom_JJA_1yr <- pdp::partial(
  object = model_sat, 
  pred.var = "tasmax_anom_JJA_1yr", 
  train = calib_dt
)
par_tasmax_anom_JJA_1yr$yhat <- exp(par_tasmax_anom_JJA_1yr$yhat) * 100
plot_tasmax_anom_JJA_1yr <- autoplot(par_tasmax_anom_JJA_1yr) +
  theme_bw() +
  ylim(c(0, 30)) +
  ggtitle("tasmax_anom_JJA_1yr [°C]") +
  ylab("Mortality rate [%]")+
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# 7 rsds_anom_MAM         
par_rsds_anom_MAM  <- pdp::partial(
  object = model_sat, 
  pred.var = "rsds_anom_MAM", 
  train = calib_dt
)

par_rsds_anom_MAM$yhat <- exp(par_rsds_anom_MAM$yhat) * 100
plot_rsds_anom_MAM <- autoplot(par_rsds_anom_MAM) +
  theme_bw() +
  ylim(c(0, 30)) +
  ggtitle(expression(rsds_anom_MAM ~ "[" * Wm^{-2} * "]")) +
  ylab("Mortality rate [%]")+
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Combine all plots
all_pdp <- ggarrange(
  plot_cwb_JJA_roll_min_5yrs,
  plot_rsds_anom_JJA_roll_mean_2yrs,
  plot_rsds_anom_MAM,
  plot_tasmin_DJF,
  plot_tasmax_anom_JJA_1yr, 
  plot_late_frosts,
  ncol = 3, nrow = 2
  
)

all_pdp <- annotate_figure(
  all_pdp,
  top = text_grob("Partial dependence plots for all predictors\n", size = 14)
)

# Display the combined plot
print(all_pdp)

ggsave("/storage/tanders/standing_deadwood_data/figures_manuscript_2/all_pdp.jpeg", plot = all_pdp, 
       width = 32, height = 16, dpi = 500, units = "cm")
ggsave("/storage/tanders/standing_deadwood_data/figures_manuscript_2/FIG3.tif", plot = all_pdp, 
       width = 32, height = 16, dpi = 500, units = "cm")

################################
### Combine the figures
################################
vip_test_gg2 <- vip_test_gg + 
  theme(plot.margin = margin(t = 10, r = 20, b = 20, l = 10))  # add space on right
all_pdp2      <- all_pdp      + theme(plot.margin = margin(t = 20))  # space above


comb_2_3_gg <- ggarrange(
  vip_test_gg2,
  all_pdp2,
  nrow = 2, ncol = 1,
  labels = c("a", "b"),
  heights = c(1, 1.5)
)

ggsave("/storage/tanders/standing_deadwood_data/figures_manuscript_2/FIG2.jpeg", plot = comb_2_3_gg, 
       width = 35, height = 40, dpi = 500, units = "cm")
ggsave("/storage/tanders/standing_deadwood_data/figures_manuscript_2/FIG2.tif", plot = comb_2_3_gg, 
       width = 35, height = 40, dpi = 500, units = "cm")


#### 8. Calibration plots ####

# Model 1
calib_dt$pred <- predict(model_sat, newdata = calib_dt, type = "response")
test_m_sat <- test_data[,c(1:4, 12)]
test_m_sat$model <- "model sat"
test_data <- na.omit(test_data)
test_m_sat <- na.omit(test_m_sat)
r_m_sat <- cor(test_m_sat$inc_deadwood_spruce, test_m_sat$pred, method = "pearson")
r2_m_sat <- r_m_sat^2
aic <- model_sat$aic

# Compute RMSE
rmse <- sqrt(mean((test_data$pred - test_data$inc_deadwood_spruce)^2))

# Normalize RMSE (using the range of observed values)
nrmse_range <- rmse / (max(test_data$inc_deadwood_spruce) - min(test_data$inc_deadwood_spruce))

# Normalize RMSE (using the mean of observed values)
nrmse_mean_test <- rmse / mean(test_data$inc_deadwood_spruce)
nrmse_mean_test

# Calculate annual mean
annual_mean <- test_m_sat[, .(mean_inc_deadwood = mean(inc_deadwood_spruce, na.rm = TRUE),
                              mean_pred = mean(pred, na.rm = TRUE)), by = year]

# Print the result
print(annual_mean)

train_m_sat <- train_data
train_m_sat$model <- "model sat"
train_data <- na.omit(train_data)
train_m_sat$tasmin_DJF <- train_m_sat$tasmin_DJF - 273.15
train_m_sat$pred <- predict(model_sat, newdata = train_m_sat, type = "response")
train_m_sat <- na.omit(train_m_sat)
r_m_sat_train <- cor(train_m_sat$inc_deadwood_spruce, train_m_sat$pred, method = "pearson")
r2_m_sat_train <- r_m_sat_train^2

# Calculate annual mean
annual_mean <- train_m_sat[, .(mean_inc_deadwood = mean(inc_deadwood_spruce, na.rm = TRUE),
                               mean_pred = mean(pred, na.rm = TRUE)), by = year]

# Print the result
print(annual_mean)

# Fit linear model
lm_fit <- lm(pred ~ inc_deadwood_spruce, data = test_m_sat)

# Extract coefficients
slope <- coef(lm_fit)[2]
intercept <- coef(lm_fit)[1]
r2_value <- summary(lm_fit)$r.squared

PI <- (slope + r2_value + dev_expl + (1-nrmse_mean_test)) / 4

# Create the calibration plot of the testing data
cal_plot_test <- ggplot(test_m_sat, aes(x=inc_deadwood_spruce*100, y=pred*100)) +
  geom_point(color = "grey50", alpha = 0.3) +
  ylim(0,75) +
  xlim(0,75) +
  xlab("Observed mortality rate (Schiefer et al. 2023) [%]") +
  ylab("Predicted mortality rate [%]") +
  geom_smooth(method = "lm", fullrange = TRUE, se = T, color = "grey20") +
  geom_abline (slope=1, linetype = "dashed", color="Red") +
  annotate(geom = "text", x = 5, y = 69, label = paste0('dev. expl. = ', round(dev_expl, 4)), colour ="black", hjust = 0) +
  annotate(geom = "text", x = 5, y = 65, label = paste0('R² = ', round(r2_m_sat, 4)), colour ="black", hjust = 0) +
  annotate(geom = "text", x = 5, y = 61, label =paste0('NRMSE = ', round(nrmse_mean_test, 4)), colour ="black", hjust = 0) +
  annotate(geom = "text", x = 5, y = 57, label =paste0('m = ', round(slope, 4)), colour ="black", hjust = 0) +
  # annotate(geom = "text", x = 5, y = 53, label =paste0('PI = ', round(PI, 4)), colour ="black", hjust = 0) +
  #  annotate(geom = "text", x = 5, y = 49, label =paste0('AIC = ', round(aic, 1)), colour ="black", hjust = 0) +
  ggtitle(paste0("Calibration plot testing data (n = ",nrow(test_data),")")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.key.size = unit(2, "line"),
        plot.subtitle = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10, vjust = +2),
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), 'cm'))
cal_plot_test

# Compute RMSE
rmse <- sqrt(mean((calib_dt$pred - calib_dt$inc_deadwood_spruce)^2))

# Normalize RMSE (using the range of observed values)
nrmse_range <- rmse / (max(calib_dt$inc_deadwood_spruce) - min(calib_dt$inc_deadwood_spruce))

# Normalize RMSE (using the mean of observed values)
nrmse_mean_train <- rmse / mean(calib_dt$inc_deadwood_spruce)
nrmse_mean_train

# Fit linear model
lm_fit <- lm(pred ~ inc_deadwood_spruce, data = train_m_sat)

# Extract coefficients
slope <- coef(lm_fit)[2]
intercept <- coef(lm_fit)[1]
r2_value <- summary(lm_fit)$r.squared

PI <- (slope + r2_value + dev_expl + (1-nrmse_mean_train)) / 4

cal_plot_train <- ggplot(train_m_sat, aes(x=inc_deadwood_spruce*100, y=pred*100)) +
  geom_point(color = "grey50", alpha = 0.3) +
  ylim(0,75) +
  xlim(0,75) +
  xlab("Observed mortality rate (Schiefer et al. 2023) [%]") +
  ylab("Predicted mortality rate [%]") +
  geom_smooth(method = "lm", fullrange = TRUE, se = T, color = "grey20") +
  geom_abline (slope=1, linetype = "dashed", color="Red") +
  #annotate("rect", xmin = 0.25, xmax = 5, ymin = 9.5, ymax = 13.6, alpha = 0.5, fill = "lightgrey") +
  #scale_color_manual(values =  c("#D55E00", "#009E73", "#56B4E9", "#CC79A7")) +
  annotate(geom = "text", x = 5, y = 69, label = paste0('dev. expl. = ', round(dev_expl, 4)), colour ="black", hjust = 0) +
  annotate(geom = "text", x = 5, y = 65, label = paste0('R² = ', round(r2_m_sat_train, 4)), colour ="black", hjust = 0) +
  annotate(geom = "text", x = 5, y = 61, label =paste0('NRMSE = ', round(nrmse_mean_train, 4)), colour ="black", hjust = 0) +
  annotate(geom = "text", x = 5, y = 57, label =paste0('m = ', round(slope, 4)), colour ="black", hjust = 0) +
  #annotate(geom = "text", x = 5, y = 53, label =paste0('PI = ', round(PI, 4)), colour ="black", hjust = 0) +
  #annotate(geom = "text", x = 5, y = 49, label =paste0('AIC = ', round(aic, 1)), colour ="black", hjust = 0) +
  ggtitle(paste0("Calibration plot training data (n = ",nrow(train_data),")")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.key.size = unit(2, "line"),
        plot.subtitle = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10, vjust = +2),
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), 'cm'))
cal_plot_train

calibration_plots <- ggarrange(cal_plot_train, cal_plot_test, labels = c("a", "b"))
ggsave("./figures_filtered/calibration_plots.png", plot = calibration_plots, width = 12, height = 6, dpi = 300, units = "in")
