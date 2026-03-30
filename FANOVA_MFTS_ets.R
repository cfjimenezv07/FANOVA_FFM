# TWA+ FTSM CFJV model

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
female_prefecture_rate_combo_smooth <- readRDS(file.path(dir.r, "female_prefecture_rate_combo_smooth.rds"))
male_prefecture_rate_combo_smooth   <- readRDS(file.path(dir.r, "male_prefecture_rate_combo_smooth.rds"))

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
# 3. CFJV2024 APPROACH: FPCA ON RESIDUALS & EXPANDING WINDOW
################################################################################
library(doMC)
registerDoMC(detectCores() - 2)

# Helper function for Eigenvalue Ratio selection (as per CFJV2024)
select_k_evr <- function(tau = 10^-2, eigenvalue) {
  k_max = length(eigenvalue)
  k_all = rep(0, k_max-1)
  for(k in 1:(k_max-1)) {
    k_all[k] = (eigenvalue[k+1]/eigenvalue[k]) * ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + 
      ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
  }
  return(which.min(k_all))
}

# Forecasting Function tailored to your 3D residual structure
# This function performs FPCA on the residuals and forecasts the scores
pref_fpca_forecast <- function(residuals_3d_slice, fh, prediction_method="ets", select_K="EVR", K=6) {
  # residuals_3d_slice is [Year x Age] for a specific prefecture
  med_polish_resi = t(residuals_3d_slice) # [Age x Year]
  med_polish_resi_cov = cov(t(med_polish_resi))
  med_polish_resi_eigen = eigen(med_polish_resi_cov)
  
  if(select_K == "Fixed") {
    retain_component = K
  } else {
    retain_component = select_k_evr(tau = 10^-2, eigenvalue = med_polish_resi_eigen$values)
  }
  
  # Basis functions and scores
  basis = as.matrix(med_polish_resi_eigen$vectors[, 1:retain_component])
  score = crossprod(med_polish_resi, basis)
  
  # Forecast PC scores
  score_fore = matrix(NA, retain_component, fh)
  for(ip in 1:retain_component) {
    # Using ets as in your reference code
    fit = ets(ts(as.vector(score[, ip])))
    score_fore[ip, ] = as.vector(forecast(fit, h = fh)$mean)
  }
  
  # Reconstruct stochastic part: [Age x Horizon]
  resi_fore = basis %*% score_fore
  return(t(resi_fore)) # Return as [Horizon x Age]
}

################################################################################
# 4. EXPANDING WINDOW EVALUATION
################################################################################

max_h <- 10
n_years <- dim(FANOVA_res_male_3D)[1] # 49
n_training_ini <- 39 # Starting training size (1975 to 2013)

# Lists to store results for each prefecture
Boot_male_pref <- vector("list", 47)
Boot_female_pref <- vector("list", 47)

for (i in 1:47) {
  cat(sprintf("Processing Prefecture %d: %s...\n", i, state_names[i]))
  
  # Extract residual slices [Year x Age]
  res_m_slice <- FANOVA_res_male_3D[, i, ]
  res_f_slice <- FANOVA_res_female_3D[, i, ]
  
  # Parallel loop over expanding windows
  out <- foreach(iwk = 1:max_h) %dopar% {
    train_end <- n_training_ini + iwk - 1
    fh_current <- n_years - train_end
    
    # 1. Forecast Residuals (Stochastic Part)
    fore_res_m <- pref_fpca_forecast(res_m_slice[1:train_end, ], fh = fh_current)
    fore_res_f <- pref_fpca_forecast(res_f_slice[1:train_end, ], fh = fh_current)
    
    # 2. Add Deterministic Part (Fixed Effects) and Back-transform
    # Deterministic part is constant over time in this model, so we slice the first 'fh_current' years
    # of the holdout period from your existing deterministic 3D arrays
    pred_rate_m <- matrix(NA, fh_current, length(ages))
    pred_rate_f <- matrix(NA, fh_current, length(ages))
    
    for(h in 1:fh_current) {
      target_yr <- train_end + h
      pred_rate_m[h, ] <- exp(FANOVA_deterministic_male_3D[target_yr, i, ] + fore_res_m[h, ])
      pred_rate_f[h, ] <- exp(FANOVA_deterministic_female_3D[target_yr, i, ] + fore_res_f[h, ])
    }
    
    list(m = pred_rate_m, f = pred_rate_f)
  }
  
  Boot_male_pref[[i]] <- lapply(out, function(x) x$m)
  Boot_female_pref[[i]] <- lapply(out, function(x) x$f)
}

################################################################################
# 5. ERROR CALCULATION (MAE & RMSE)
################################################################################

# Initialize lists for all metrics per horizon
err_m_mae  <- err_f_mae  <- lapply(1:max_h, function(x) numeric())
err_m_rmse <- err_f_rmse <- lapply(1:max_h, function(x) numeric())

for (i in 1:47) {
  for (iwk in 1:max_h) {
    train_end <- n_training_ini + iwk - 1
    fh_current <- n_years - train_end
    
    for (h in 1:fh_current) {
      target_yr <- train_end + h
      
      # 1. Get Truth and Forecast
      actual_m <- truth_male_3d[, i, target_yr]
      actual_f <- truth_female_3d[, i, target_yr]
      
      pred_m <- Boot_male_pref[[i]][[iwk]][h, ]
      pred_f <- Boot_female_pref[[i]][[iwk]][h, ]
      
      # 2. Compute MAE
      err_m_mae[[h]] <- c(err_m_mae[[h]], mean(abs(pred_m - actual_m)))
      err_f_mae[[h]] <- c(err_f_mae[[h]], mean(abs(pred_f - actual_f)))
      
      # 3. Compute RMSE
      err_m_rmse[[h]] <- c(err_m_rmse[[h]], sqrt(mean((pred_m - actual_m)^2)))
      err_f_rmse[[h]] <- c(err_f_rmse[[h]], sqrt(mean((pred_f - actual_f)^2)))
    }
  }
}

# Final Summary Table
final_summary <- data.frame(
  Horizon     = 1:max_h,
  Male_MAE    = sapply(err_m_mae, mean),
  Male_RMSE   = sapply(err_m_rmse, mean),
  Female_MAE  = sapply(err_f_mae, mean),
  Female_RMSE = sapply(err_f_rmse, mean)
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

