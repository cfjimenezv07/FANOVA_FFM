# 1. Two-Way Anova + FFM (Omitting One-Way ANOVA)

# FANOVA_FFM forecasting evaluation
setwd("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Rcodes")

# dir.r for results
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Results"
# dir.p
dir.p <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Plots"

################################################################################
# 1. SETUP & DATA PREPARATION
################################################################################
source("1_load_packages.R")
source("FFM.R") # Contains HDFTS_factor_decomp

# Prefecture Names in Order
state_names <- c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
                 "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
                 "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", 
                 "Shizuoka", "Aichi", "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", 
                 "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
                 "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", 
                 "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

n_pref <- length(state_names) # 47
ages   <- 0:95
years  <- 1975:2023
n_years <- length(years)      # 49

# Load Data

# female_smooth / male_smooth: [96 x 2303] (Ages x (Years*Prefectures))
female_prefecture_rate_combo_smooth   <- readRDS(file.path(dir.r, "female_prefecture_rate_combo_smooth.rds"))
male_prefecture_rate_combo_smooth     <- readRDS(file.path(dir.r, "male_prefecture_rate_combo_smooth.rds"))
# truth_3d: [96 x 47 x 49] (Age x Year x Prefecture)
truth_female_3d = readRDS(file = file.path(dir.r, "female_prefecture_rate_array.rds"))
truth_male_3d = readRDS(file = file.path(dir.r, "male_prefecture_rate_array.rds"))

################################################################################
# 2. TWO-WAY ANOVA DECOMPOSITION
################################################################################

# Compute Functional Means (TWA)
FANOVA_means <- FANOVA(
  data_pop1 = log(male_prefecture_rate_combo_smooth),
  data_pop2 = log(female_prefecture_rate_combo_smooth),
  year = years, age = ages, n_prefectures = 47, n_populations = 2
)

# Extract Residuals directly from TWA
FANOVA_residuals_raw <- Two_way_Residuals_means(
  data_pop1 = log(male_prefecture_rate_combo_smooth), 
  data_pop2 = log(female_prefecture_rate_combo_smooth), 
  year = years, age = ages, n_prefectures = 47, n_populations = 2
)

# Reshape TWA residuals to 3D for Factor Model [Year x Prefecture x Age]
FANOVA_res_female_3D <- array(FANOVA_residuals_raw$residuals2_mean, dim = c(length(years), 47, 96))
FANOVA_res_male_3D   <- array(FANOVA_residuals_raw$residuals1_mean, dim = c(length(years), 47, 96))

# Construct Deterministic Part (TWA Fixed Effects Only)
# We expand the TWA components to match the [Year x Prefecture x Age] structure
FANOVA_deterministic_male_3D   <- array(NA, dim = c(49, 47, 96))
FANOVA_deterministic_female_3D <- array(NA, dim = c(49, 47, 96))

for(i in 1:47) {
  # Deterministic = Grand Effect + Region Effect + Cohort (Population) Effect
  det_val_m <- FANOVA_means$FGE_mean + FANOVA_means$FRE_mean[i,] + FANOVA_means$FCE_mean[1,]
  det_val_f <- FANOVA_means$FGE_mean + FANOVA_means$FRE_mean[i,] + FANOVA_means$FCE_mean[2,]
  
  for(t in 1:49) {
    FANOVA_deterministic_male_3D[t, i, ]   <- det_val_m
    FANOVA_deterministic_female_3D[t, i, ] <- det_val_f
  }
}

################################################################################
# 3. FUNCTIONAL FACTOR MODEL (Applied to TWA Residuals)
################################################################################

HDFTS_m <- HDFTS_factor_decomp(FANOVA_res_male_3D)
HDFTS_f <- HDFTS_factor_decomp(FANOVA_res_female_3D)

factors_m_all  <- HDFTS_m$factors
loadings_m_all <- HDFTS_m$factor_loading # [K, Age, Prefecture]

factors_f_all  <- HDFTS_f$factors
loadings_f_all <- HDFTS_f$factor_loading

################################################################################
# 4. EXPANDING WINDOW (EW) EVALUATION
################################################################################

err_m_mae_list  <- err_f_mae_list  <- lapply(1:10, function(x) numeric())
err_m_rmse_list <- err_f_rmse_list <- lapply(1:10, function(x) numeric())

for (ew in 1:10) {
  train_end <- 38 + ew 
  fh_max    <- 49 - train_end
  
  cat(sprintf("EW %d: Forecasting %d horizons...\n", ew, fh_max))
  
  # Forecast Factors
  f_fore_m_raw <- sapply(1:ncol(factors_m_all), function(k) {
    forecast(ets(ts(factors_m_all[1:train_end, k])), h = fh_max)$mean
  })
  f_fore_f_raw <- sapply(1:ncol(factors_f_all), function(k) {
    forecast(ets(ts(factors_f_all[1:train_end, k])), h = fh_max)$mean
  })
  
  f_fore_m <- matrix(f_fore_m_raw, nrow = fh_max, ncol = ncol(factors_m_all))
  f_fore_f <- matrix(f_fore_f_raw, nrow = fh_max, ncol = ncol(factors_f_all))
  
  # Reconstruct and Evaluate
  for (h in 1:fh_max) {
    target_idx <- train_end + h
    p_mae_m <- p_mae_f <- p_rmse_m <- p_rmse_f <- numeric(47)
    
    f_h_m <- f_fore_m[h, , drop = FALSE]
    f_h_f <- f_fore_f[h, , drop = FALSE]
    
    for (p in 1:47) {
      # 1. Reconstruct Stochastic Part from Factors
      res_hat_m <- as.numeric(f_h_m %*% loadings_m_all[, , p])
      res_hat_f <- as.numeric(f_h_f %*% loadings_f_all[, , p])
      
      # 2. Get Deterministic Part (TWA Fixed Effects)
      det_m <- FANOVA_deterministic_male_3D[target_idx, p, ]
      det_f <- FANOVA_deterministic_female_3D[target_idx, p, ]
      
      # 3. Back-transform Prediction 
      pred_m <- exp(det_m + res_hat_m)
      pred_f <- exp(det_f + res_hat_f)
      
      # 4. Get Truth
      actual_m <- truth_male_3d[, p, target_idx]
      actual_f <- truth_female_3d[, p, target_idx]
      
      # 5. Compute Metrics
      p_mae_m[p]  <- ftsa:::mae(pred_m, actual_m)
      p_rmse_m[p] <- ftsa:::rmse(pred_m, actual_m)
      p_mae_f[p]  <- ftsa:::mae(pred_f, actual_f)
      p_rmse_f[p] <- ftsa:::rmse(pred_f, actual_f)
    }
    
    err_m_mae_list[[h]]  <- c(err_m_mae_list[[h]],  mean(p_mae_m))
    err_m_rmse_list[[h]] <- c(err_m_rmse_list[[h]], mean(p_rmse_m))
    err_f_mae_list[[h]]  <- c(err_f_mae_list[[h]],  mean(p_mae_f))
    err_f_rmse_list[[h]] <- c(err_f_rmse_list[[h]], mean(p_rmse_f))
  }
}

# Final Summary Table
final_summary <- data.frame(
  fh = 1:10, 
  Male_MAE    = sapply(err_m_mae_list,  mean), 
  Male_RMSE   = sapply(err_m_rmse_list, mean),
  Female_MAE  = sapply(err_f_mae_list,  mean),
  Female_RMSE = sapply(err_f_rmse_list, mean)
)

print(final_summary)
# 1. Calculate the mean and median for all metric columns (cols 2 to 5)
summary_means   <- colMeans(final_summary[, -1], na.rm = TRUE)
summary_medians <- apply(final_summary[, -1], 2, median, na.rm = TRUE)

# 2. Convert fh to character so we can add the "Mean"/"Median" labels
final_summary$fh <- as.character(final_summary$fh)

# 3. Create the new rows as data frames to match the structure
mean_row   <- data.frame(fh = "Mean",   t(summary_means))
median_row <- data.frame(fh = "Median", t(summary_medians))

# 4. Use rbind to append them to the bottom
final_summary <- rbind(final_summary, mean_row, median_row)

# Optional: ensure metrics are numeric if needed for further calculations
final_summary[, 2:5] <- lapply(final_summary[, 2:5], as.numeric)

print(final_summary)
saveRDS(final_summary, file.path(dir.r, "TWA_FFM_final_summary_ets.rds"))


#ARIMA


# Initialize lists for MAE and RMSE
err_m_mae_list <- err_f_mae_list <- lapply(1:10, function(x) numeric())
err_m_rmse_list <- err_f_rmse_list <- lapply(1:10, function(x) numeric())

for (ew in 1:10) {
  train_end <- 38 + ew 
  fh_max    <- 49 - train_end
  
  cat(sprintf("EW %d: Forecasting %d horizons using ARIMA...\n", ew, fh_max))
  
  # --- Forecast Factors using auto.arima ---
  f_fore_m_raw <- sapply(1:ncol(factors_m_all), function(k) {
    # Replaced ets with auto.arima
    forecast(auto.arima(ts(factors_m_all[1:train_end, k])), h = fh_max)$mean
  })
  f_fore_f_raw <- sapply(1:ncol(factors_f_all), function(k) {
    # Replaced ets with auto.arima
    forecast(auto.arima(ts(factors_f_all[1:train_end, k])), h = fh_max)$mean
  })
  
  # Ensure dimensions are [h, K]
  f_fore_m <- matrix(f_fore_m_raw, nrow = fh_max, ncol = ncol(factors_m_all))
  f_fore_f <- matrix(f_fore_f_raw, nrow = fh_max, ncol = ncol(factors_f_all))
  
  # --- Reconstruct and Evaluate ---
  for (h in 1:fh_max) {
    target_idx <- train_end + h
    p_mae_m <- p_mae_f <- p_rmse_m <- p_rmse_f <- numeric(47)
    
    f_h_m <- f_fore_m[h, , drop = FALSE]
    f_h_f <- f_fore_f[h, , drop = FALSE]
    
    for (p in 1:47) {
      # 1. Reconstruct Stochastic Part
      res_hat_m <- as.numeric(f_h_m %*% loadings_m_all[, , p])
      res_hat_f <- as.numeric(f_h_f %*% loadings_f_all[, , p])
      
      # 2. Get Deterministic Part
      det_m <- FANOVA_deterministic_male_3D[target_idx, p, ]
      det_f <- FANOVA_deterministic_female_3D[target_idx, p, ]
      
      # 3. Back-transform Prediction 
      pred_m <- exp(det_m + res_hat_m)
      pred_f <- exp(det_f + res_hat_f)
      
      # 4. Get Truth
      actual_m <- truth_male_3d[, p, target_idx]
      actual_f <- truth_female_3d[, p, target_idx]
      
      # 5. Compute Metrics
      p_mae_m[p]  <- ftsa:::mae(pred_m, actual_m)
      p_rmse_m[p] <- ftsa:::rmse(pred_m, actual_m)
      p_mae_f[p]  <- ftsa:::mae(pred_f, actual_f)
      p_rmse_f[p] <- ftsa:::rmse(pred_f, actual_f)
    }
    
    # Store results for horizon 'h'
    err_m_mae_list[[h]]  <- c(err_m_mae_list[[h]],  mean(p_mae_m))
    err_m_rmse_list[[h]] <- c(err_m_rmse_list[[h]], mean(p_rmse_m))
    err_f_mae_list[[h]]  <- c(err_f_mae_list[[h]],  mean(p_mae_f))
    err_f_rmse_list[[h]] <- c(err_f_rmse_list[[h]], mean(p_rmse_f))
  }
}

# --- Final Summary Table Construction ---
final_summary <- data.frame(
  fh = 1:10, 
  Male_MAE    = sapply(err_m_mae_list,  mean), 
  Male_RMSE   = sapply(err_m_rmse_list, mean),
  Female_MAE  = sapply(err_f_mae_list,  mean),
  Female_RMSE = sapply(err_f_rmse_list, mean)
)

# 1. Calculate Mean and Median
summary_means   <- colMeans(final_summary[, -1], na.rm = TRUE)
summary_medians <- apply(final_summary[, -1], 2, median, na.rm = TRUE)

# 2. Add labels
final_summary$fh <- as.character(final_summary$fh)
mean_row   <- data.frame(fh = "Mean",   t(summary_means))
median_row <- data.frame(fh = "Median", t(summary_medians))

# 3. Combine
final_summary <- rbind(final_summary, mean_row, median_row)

# 4. Convert metric columns back to numeric
final_summary[, 2:5] <- lapply(final_summary[, 2:5], as.numeric)

print(final_summary)
# saveRDS(final_summary, file.path(dir.r, "TWA_FFM_final_summary_arima.rds"))

