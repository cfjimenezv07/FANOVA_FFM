################################################################################
# TWA + Local FPCA (JCGS/MFTS Method) - Full Interval Evaluation Suite
################################################################################

# 1. SETUP & DATA PREPARATION
################################################################################
setwd("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Rcodes")
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Results"

source("1_load_packages.R")
source("choice_K.R")       
source("interval_score.R")  

library(forecast)
library(doMC)
registerDoMC(detectCores() - 2)

# Load Data
truth_female_3d  <- readRDS(file.path(dir.r, "female_prefecture_rate_array.rds"))
truth_male_3d    <- readRDS(file.path(dir.r, "male_prefecture_rate_array.rds"))
TWA_res_female   <- readRDS(file.path(dir.r, "TWA_res_female_3d.rds"))
TWA_res_male     <- readRDS(file.path(dir.r, "TWA_res_male_3d.rds"))
TWA_Fixed_male   <- readRDS(file.path(dir.r, "TWA_Fixed_male.rds"))
TWA_Fixed_female <- readRDS(file.path(dir.r, "TWA_Fixed_female.rds"))

n_pref <- 47; ages <- 0:95; years <- 1975:2023; n_years <- length(years)

################################################################################
# 2. CORE JCGS (MFTS) FUNCTIONS
################################################################################
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
# Local FPCA Forecasting Engine (JCGS / MFTS Method)
# This processes each prefecture as an independent functional time series
JCGS_MFTS_forecast_fun <- function(data_set_array, fh, fmethod, Fixed_TWA) {
  n_age = dim(data_set_array)[1]
  n_state = dim(data_set_array)[2]
  n_year = dim(data_set_array)[3]
  
  forecast_matrix = matrix(NA, n_age, n_state)
  
  for(i in 1:n_state) {
    # Local slice for one prefecture
    slice = data_set_array[, i, ]
    med_polish_resi_cov = cov(t(slice))
    eigen_decomp = eigen(med_polish_resi_cov)
    
    # Selection of components K per prefecture
    retain_K = select_k_evr(tau = 10^-2, eigenvalue = eigen_decomp$values)
    
    basis = as.matrix(eigen_decomp$vectors[, 1:retain_K])
    scores = crossprod(slice, basis)
    
    score_fore = numeric(retain_K)
    for(k in 1:retain_K) {
      if(fmethod == "ets") {
        score_fore[k] = forecast(ets(scores[, k]), h = fh)$mean[fh]
      } else {
        score_fore[k] = forecast(auto.arima(scores[, k]), h = fh)$mean[fh]
      }
    }
    # Stochastic Reconstruction + Local Fixed Effects
    forecast_matrix[, i] = (basis %*% score_fore) + Fixed_TWA[, i]
  }
  return(forecast_matrix)
}

# Calibration Helper
tune_para_select <- function(resi, tune_para, nominal, resi_sd) {
  n_age = length(resi_sd)
  ind = ifelse(abs(resi) <= tune_para * resi_sd, 1, 0)
  return(abs(sum(ind)/(n_age * ncol(resi)) - nominal))
}

# Interval Evaluation Function (MFTS Logic)
FM_TWA_JCGS_fun_int_eval <- function(raw_data, smooth_data, fh, fore_method, alpha, Fixed_TWA) {
  n_age = dim(raw_data)[1]; n_state = dim(raw_data)[2]; n_year = dim(raw_data)[3]
  
  # 1. Validation for calibration
  smooth_mat = array(NA, dim = c(n_age, n_state, (12 - fh)))
  for(iw in 1:(12 - fh)) {
    smooth_mat[,,iw] <- JCGS_MFTS_forecast_fun(smooth_data[,,1:((n_year - 22) + iw)], fh, fore_method, Fixed_TWA)
  }
  
  # Forecast Errors in original scale
  err_mat <- array(raw_data[,,(n_year - (21 - fh)):(n_year - 10)], dim = c(n_age, n_state, (12 - fh))) - exp(smooth_mat)
  
  MFTS_quantile = MFTS_sd = matrix(NA, n_age, n_state)
  tune_para_value = vector("numeric", n_state)
  
  for(iw in 1:n_state) {
    MFTS_quantile[,iw] = apply(abs(err_mat[,iw,]), 1, quantile, (1 - alpha))
    MFTS_sd[,iw] = apply(err_mat[,iw,], 1, sd)
    
    # Optimization for parameter 'c' (tune_para)
    obj_BFGS = optim(par = 1, fn = tune_para_select, method = "BFGS", resi = err_mat[,iw,], nominal = 1 - alpha, resi_sd = MFTS_sd[,iw])
    obj_NM   = optim(par = 1, fn = tune_para_select, method = "Nelder-Mead", resi = err_mat[,iw,], nominal = 1 - alpha, resi_sd = MFTS_sd[,iw])
    obj_Brent= optim(par = 1, fn = tune_para_select, method = "Brent", resi = err_mat[,iw,], nominal = 1 - alpha, lower = 0, upper = 10, resi_sd = MFTS_sd[,iw])
    
    tune_para_value[iw] = c(obj_BFGS$par, obj_NM$par, obj_Brent$par)[which.min(c(obj_BFGS$value, obj_NM$value, obj_Brent$value))]
  }
  
  # 2. Test Set Forecasting
  smooth_mat_fore = array(NA, dim = c(n_age, n_state, (11 - fh)))
  for(iw in 1:(11 - fh)) {
    smooth_mat_fore[,,iw] <- JCGS_MFTS_forecast_fun(smooth_data[,,1:((n_year - 11) + iw)], fh, fore_method, Fixed_TWA)
  }
  
  point_forecast <- exp(smooth_mat_fore)
  int_score_sd = int_score_quantile = matrix(NA, n_state, 3)
  
  for(iw in 1:n_state) {
    # SD-based Intervals
    lb_sd = point_forecast[,iw,] - tune_para_value[iw] * MFTS_sd[,iw]
    ub_sd = point_forecast[,iw,] + tune_para_value[iw] * MFTS_sd[,iw]
    
    # Quantile-based Intervals
    lb_quant = point_forecast[,iw,] - MFTS_quantile[,iw]
    ub_quant = point_forecast[,iw,] + MFTS_quantile[,iw]
    
    actuals = raw_data[,iw,(n_year - (10 - fh)):n_year]
    int_score_sd[iw,] = interval_score(holdout = actuals, lb = lb_sd, ub = ub_sd, alpha = alpha)
    int_score_quantile[iw,] = interval_score(holdout = actuals, lb = lb_quant, ub = ub_quant, alpha = alpha)
  }
  return(list(int_score_sd = int_score_sd, int_score_quantile = int_score_quantile))
}

################################################################################
# 3. AUTOMATED EXECUTION HANDLER
################################################################################

run_scenario <- function(gender, raw, smooth, fixed, alpha_val, method) {
  cat(sprintf("\nRunning MFTS: %s | %s | Alpha: %s\n", gender, method, alpha_val))
  
  res_sd = res_quant = array(NA, dim = c(47, 3, 10), dimnames = list(1:47, c("ECP", "CPD", "score"), 1:10))
  
  for(h in 1:10) {
    cat("H:", h, " ")
    dum = FM_TWA_JCGS_fun_int_eval(raw, smooth, h, method, alpha_val, fixed)
    res_sd[,,h] = dum$int_score_sd
    res_quant[,,h] = dum$int_score_quantile
  }
  
  fn <- sprintf("JCGS_MFTS_%s_%s_%s", tolower(gender), method, alpha_val)
  saveRDS(res_sd, file.path(dir.r, paste0(fn, "_sd.rds")))
  saveRDS(res_quant, file.path(dir.r, paste0(fn, "_quant.rds")))
  saveRDS(t(apply(res_sd, c(2, 3), mean)), file.path(dir.r, paste0(fn, "_sd_mean.rds")))
  saveRDS(t(apply(res_quant, c(2, 3), mean)), file.path(dir.r, paste0(fn, "_quant_mean.rds")))
}

################################################################################
# 4. FULL EXECUTION (All combinations)
################################################################################

# FEMALE
run_scenario("Female", truth_female_3d, TWA_res_female, TWA_Fixed_female, 0.20, "arima")
run_scenario("Female", truth_female_3d, TWA_res_female, TWA_Fixed_female, 0.20, "ets")
run_scenario("Female", truth_female_3d, TWA_res_female, TWA_Fixed_female, 0.05, "arima")
run_scenario("Female", truth_female_3d, TWA_res_female, TWA_Fixed_female, 0.05, "ets")

# MALE
run_scenario("Male", truth_male_3d, TWA_res_male, TWA_Fixed_male, 0.20, "arima")
run_scenario("Male", truth_male_3d, TWA_res_male, TWA_Fixed_male, 0.20, "ets")
run_scenario("Male", truth_male_3d, TWA_res_male, TWA_Fixed_male, 0.05, "arima")
run_scenario("Male", truth_male_3d, TWA_res_male, TWA_Fixed_male, 0.05, "ets")

cat("\n*** JCGS MFTS EVALUATION COMPLETE ***\n")


################################################################################
# FULL SUMMARY: JCGS MFTS 0.05 RESULTS WITH MEAN & MEDIAN ROWS
################################################################################

# 1. HELPER FUNCTION
################################################################################
# This function calculates Mean/Median and appends them as the 11th and 12th rows
append_stats <- function(df) {
  # Calculate Mean and Median down the columns (ECP, CPD, score)
  stats_mean   <- apply(df, 2, mean)
  stats_median <- apply(df, 2, median)
  
  # Combine with the original 10-horizon table
  final_tab <- rbind(df, Mean = stats_mean, Median = stats_median)
  return(final_tab)
}

# 2. PROCESSING ALL 8 SCENARIOS (Alpha = 0.05)
################################################################################

# --- FEMALE ARIMA ---
mean_sd_female_arima_0.05    <- append_stats(readRDS(file.path(dir.r, "JCGS_MFTS_female_arima_0.05_sd_mean.rds")))
mean_quant_female_arima_0.05 <- append_stats(readRDS(file.path(dir.r, "JCGS_MFTS_female_arima_0.05_quant_mean.rds")))

# --- FEMALE ETS ---
mean_sd_female_ets_0.05      <- append_stats(readRDS(file.path(dir.r, "JCGS_MFTS_female_ets_0.05_sd_mean.rds")))
mean_quant_female_ets_0.05    <- append_stats(readRDS(file.path(dir.r, "JCGS_MFTS_female_ets_0.05_quant_mean.rds")))

# --- MALE ARIMA ---
mean_sd_male_arima_0.05      <- append_stats(readRDS(file.path(dir.r, "JCGS_MFTS_male_arima_0.05_sd_mean.rds")))
mean_quant_male_arima_0.05    <- append_stats(readRDS(file.path(dir.r, "JCGS_MFTS_male_arima_0.05_quant_mean.rds")))

# --- MALE ETS ---
mean_sd_male_ets_0.05        <- append_stats(readRDS(file.path(dir.r, "JCGS_MFTS_male_ets_0.05_sd_mean.rds")))
mean_quant_male_ets_0.05      <- append_stats(readRDS(file.path(dir.r, "JCGS_MFTS_male_ets_0.05_quant_mean.rds")))


# 3. PRINT ALL RESULTS
################################################################################

cat("\n========================================================\n")
cat("FEMALE - ALPHA 0.05 - INTERVAL EVALUATION (AVERAGED)\n")
cat("========================================================\n")

cat("\n[1] ARIMA - SD Intervals:\n")
print(round(mean_sd_female_arima_0.05, 7))

cat("\n[2] ARIMA - Quantile Intervals:\n")
print(round(mean_quant_female_arima_0.05, 7))

cat("\n[3] ETS - SD Intervals:\n")
print(round(mean_sd_female_ets_0.05, 7))

cat("\n[4] ETS - Quantile Intervals:\n")
print(round(mean_quant_female_ets_0.05, 7))

cat("\n\n========================================================\n")
cat("MALE - ALPHA 0.05 - INTERVAL EVALUATION (AVERAGED)\n")
cat("========================================================\n")

cat("\n[5] ARIMA - SD Intervals:\n")
print(round(mean_sd_male_arima_0.05, 7))

cat("\n[6] ARIMA - Quantile Intervals:\n")
print(round(mean_quant_male_arima_0.05, 7))

cat("\n[7] ETS - SD Intervals:\n")
print(round(mean_sd_male_ets_0.05, 7))

cat("\n[8] ETS - Quantile Intervals:\n")
print(round(mean_quant_male_ets_0.05, 7))

cat("\n--- ALL TABLES PRINTED ---\n")
