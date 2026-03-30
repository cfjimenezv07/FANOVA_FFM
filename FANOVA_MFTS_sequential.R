################################################################################
# COMPLETE JCGS/MFTS SEQUENTIAL EVALUATION 
################################################################################

# 1. SETUP & PATHS
################################################################################
setwd("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Rcodes")
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Results"

source("1_load_packages.R")
source("interval_score.R") # Ensure this file contains your interval_score() function

library(quantreg)
library(forecast)
library(doMC)
registerDoMC(detectCores() - 2)

# 2. AUXILIARY FUNCTIONS (JCGS / MFTS LOGIC)
################################################################################

# EVR Criterion for selecting K components per prefecture
select_k_evr <- function(tau = 10^-2, eigenvalue) {
  k_max = length(eigenvalue)
  k_all = rep(0, k_max-1)
  for(k in 1:(k_max-1)) {
    k_all[k] = (eigenvalue[k+1]/eigenvalue[k]) * ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + 
      ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
  }
  return(which.min(k_all))
}

# Local FPCA Forecasting Engine (The JCGS approach)
FM_TWA_JCGS_fun <- function(data_set_array, fh, fmethod, Fixed_TWA) {
  n_age = dim(data_set_array)[1]
  n_state = dim(data_set_array)[2]
  n_year = dim(data_set_array)[3]
  
  JCGS_forecast = matrix(NA, n_age, n_state)
  
  for(i in 1:n_state) {
    # Extract local slice for one prefecture
    slice = data_set_array[, i, ]
    med_polish_resi_cov = cov(t(slice))
    eigen_decomp = eigen(med_polish_resi_cov)
    
    # Select K locally
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
    # Reconstruction + Prefecture Fixed Effect
    JCGS_forecast[, i] = (basis %*% score_fore) + Fixed_TWA[, i]
  }
  return(JCGS_forecast)
}

# Sequential Evaluation Function (Quantile Regression Logic)
FM_TWA_JCGS_fun_int_eval_sequential <- function(raw_data, smooth_data, fh, fore_method, alpha, Fixed_TWA) {
  n_age = dim(raw_data)[1]; n_state = dim(raw_data)[2]; n_year = dim(raw_data)[3]
  n_validation = n_year - 11
  
  int_score_array = array(NA, dim = c((11 - fh), 3, n_state), 
                          dimnames = list(1:(11 - fh), c("ECP", "CPD", "score"), 1:n_state))
  
  for(iwk in 1:(11 - fh)) {
    validation_set = 2:(n_validation + iwk - fh)
    fore_validation = array(NA, dim = c(n_age, n_state, length(validation_set)))
    
    for(ij in validation_set) {
      fore_validation[,,(ij - 1)] = FM_TWA_JCGS_fun(smooth_data[,,1:ij], fh, fore_method, Fixed_TWA)
    }
    
    # Absolute residuals in original scale
    resi = abs(raw_data[,,(2 + fh):((n_validation + iwk))] - exp(fore_validation))
    q_alpha = matrix(NA, n_age, n_state)
    
    for(iw in 1:n_state) {
      for(ij in 1:n_age) {
        # RQ Fit
        rq_fit <- rq(resi[ij,iw,] ~ 1, tau = (1 - alpha))
        # AR Order Selection
        AR_p <- max(ar(resi[ij,iw,])$order, 1)
        lastlags <- tail(x = resi[ij,iw,], AR_p)
        newdf <- as.data.frame(t(lastlags))
        names(newdf) <- paste0("lag", 1:AR_p)
        q_alpha[ij,iw] = predict(rq_fit, newdf)
      }
    }
    
    # Final Forecast for evaluation
    point_forecast <- exp(FM_TWA_JCGS_fun(smooth_data[,,1:(n_validation + iwk)], fh, fore_method, Fixed_TWA))
    
    for(iw in 1:n_state) {
      int_score_array[iwk,,iw] <- interval_score(holdout = raw_data[,iw,(n_validation + iwk + fh)], 
                                                 lb = (point_forecast - q_alpha)[,iw],
                                                 ub = (point_forecast + q_alpha)[,iw], 
                                                 alpha = alpha)
    }
  }
  return(int_score_array)
}

# 3. DATA LOADING
################################################################################
truth_female_3d  <- readRDS(file.path(dir.r, "female_prefecture_rate_array.rds"))
truth_male_3d    <- readRDS(file.path(dir.r, "male_prefecture_rate_array.rds"))
TWA_res_female   <- readRDS(file.path(dir.r, "TWA_res_female_3d.rds"))
TWA_res_male     <- readRDS(file.path(dir.r, "TWA_res_male_3d.rds"))
TWA_Fixed_male   <- readRDS(file.path(dir.r, "TWA_Fixed_male.rds"))
TWA_Fixed_female <- readRDS(file.path(dir.r, "TWA_Fixed_female.rds"))

# 4. EXECUTION HANDLER
################################################################################
run_and_save_seq <- function(gender, raw, smooth, fixed, alpha, method) {
  cat(sprintf("\n--- Start: %s | %s | Alpha: %s ---\n", gender, method, alpha))
  full_results_list = list()
  
  for(h in 1:10) {
    cat("Horizon:", h, " ")
    full_results_list[[h]] = FM_TWA_JCGS_fun_int_eval_sequential(raw, smooth, h, method, alpha, fixed)
  }
  
  # File Names
  prefix <- sprintf("JCGS_seq_%s_%s_%s", tolower(gender), method, alpha)
  
  # Save Full List
  saveRDS(full_results_list, file.path(dir.r, paste0(prefix, "_list.rds")))
  
  # Save Summary Matrix (10 x 3)
  summary_mat <- matrix(NA, 10, 3, dimnames = list(1:10, c("ECP", "CPD", "score")))
  for(h in 1:10) { summary_mat[h,] <- apply(full_results_list[[h]], 2, mean) }
  saveRDS(summary_mat, file.path(dir.r, paste0(prefix, "_mat.rds")))
  
  # Save Prefecture-specific array if alpha = 0.05
  if(alpha == 0.05) {
    pref_arr <- array(NA, dim = c(47, 3, 10))
    for(h in 1:10) { pref_arr[,,h] <- t(apply(full_results_list[[h]], c(2, 3), mean)) }
    saveRDS(pref_arr, file.path(dir.r, paste0(prefix, "_val.rds")))
  }
}

# 5. RUN ALL SCENARIOS
################################################################################

# FEMALE
run_and_save_seq("Female", truth_female_3d, TWA_res_female, TWA_Fixed_female, 0.20, "arima")
run_and_save_seq("Female", truth_female_3d, TWA_res_female, TWA_Fixed_female, 0.20, "ets")
run_and_save_seq("Female", truth_female_3d, TWA_res_female, TWA_Fixed_female, 0.05, "arima")
run_and_save_seq("Female", truth_female_3d, TWA_res_female, TWA_Fixed_female, 0.05, "ets")

# MALE
run_and_save_seq("Male", truth_male_3d, TWA_res_male, TWA_Fixed_male, 0.20, "arima")
run_and_save_seq("Male", truth_male_3d, TWA_res_male, TWA_Fixed_male, 0.20, "ets")
run_and_save_seq("Male", truth_male_3d, TWA_res_male, TWA_Fixed_male, 0.05, "arima")
run_and_save_seq("Male", truth_male_3d, TWA_res_male, TWA_Fixed_male, 0.05, "ets")

cat("\n*** ALL SEQUENTIAL PROCESSING FINISHED ***\n")



# 1. HELPER FUNCTION (Same as your provided logic)
################################################################################
append_stats <- function(df) {
  stats_mean   <- colMeans(df)
  stats_median <- apply(df, 2, median)
  
  final_tab <- rbind(df, Mean = stats_mean, Median = stats_median)
  return(final_tab)
}

# 2. PROCESSING SEQUENTIAL SCENARIOS (Alpha = 0.05)
################################################################################

# --- FEMALE ---
seq_female_arima_0.05 <- append_stats(readRDS(file.path(dir.r, "JCGS_seq_female_arima_0.05_mat.rds")))
seq_female_ets_0.05   <- append_stats(readRDS(file.path(dir.r, "JCGS_seq_female_ets_0.05_mat.rds")))

# --- MALE ---
seq_male_arima_0.05   <- append_stats(readRDS(file.path(dir.r, "JCGS_seq_male_arima_0.05_mat.rds")))
seq_male_ets_0.05     <- append_stats(readRDS(file.path(dir.r, "JCGS_seq_male_ets_0.05_mat.rds")))


# 3. PRINT SEQUENTIAL RESULTS
################################################################################

cat("\n========================================================\n")
cat("SEQUENTIAL - ALPHA 0.05 - INTERVAL EVALUATION (AVERAGED)\n")
cat("========================================================\n")

cat("\n[1] FEMALE - ARIMA Sequential:\n")
print(round(seq_female_arima_0.05, 7))

cat("\n[2] FEMALE - ETS Sequential:\n")
print(round(seq_female_ets_0.05, 7))

cat("\n--------------------------------------------------------\n")

cat("\n[3] MALE - ARIMA Sequential:\n")
print(round(seq_male_arima_0.05, 7))

cat("\n[4] MALE - ETS Sequential:\n")
print(round(seq_male_ets_0.05, 7))

cat("\n--- ALL SEQUENTIAL TABLES PRINTED ---\n")
