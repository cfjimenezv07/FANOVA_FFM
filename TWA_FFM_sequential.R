# Intervals TWA+FFM sequential

# FANOVA_FFM forecasting evaluation
setwd("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Rcodes")

#dir.r for results
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Results"
# dir.p
dir.p <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Plots"

################################################################################
# 1. SETUP & DATA PREPARATION
################################################################################
source("1_load_packages.R")
source("FFM.R")
source("choice_K.R")
source("interval_score.R")


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


# Compute Functional Means (TWA)
FANOVA_means <- FANOVA(
  data_pop1 = log(male_prefecture_rate_combo_smooth),
  data_pop2 = log(female_prefecture_rate_combo_smooth),
  year = years, age = ages, n_prefectures = 47, n_populations = 2
)

FANOVA_residuals_raw <- Two_way_Residuals_means(
  data_pop1 = log(male_prefecture_rate_combo_smooth), 
  data_pop2 = log(female_prefecture_rate_combo_smooth), 
  year = years, age = ages, n_prefectures = 47, n_populations = 2
)

FANOVA_res_female_3d <- array(FANOVA_residuals_raw$residuals2_mean, dim = c(length(years), 47, 96))
FANOVA_res_male_3d   <- array(FANOVA_residuals_raw$residuals1_mean, dim = c(length(years), 47, 96))


FANOVA_matrix_female <- matrix(aperm(FANOVA_res_female_3d, c(3, 1, 2)), nrow = 96)
FANOVA_matrix_male   <- matrix(aperm(FANOVA_res_male_3d, c(3, 1, 2)), nrow = 96)

# Load Data

# Reading using the variable dir.r
TWA_res_female <- readRDS(file.path(dir.r, "TWA_res_female_3d.rds"))
TWA_res_male   <- readRDS(file.path(dir.r, "TWA_res_male_3d.rds"))

# Read the TWA_Fixed from the TWA decomposition
TWA_Fixed_male   <- readRDS(file.path(dir.r, "TWA_Fixed_male.rds"))
TWA_Fixed_female <- readRDS(file.path(dir.r, "TWA_Fixed_female.rds"))

# truth_3d: [96 x 47 x 49] (Age x Year x Prefecture)
truth_female_3d = readRDS(file = file.path(dir.r, "female_prefecture_rate_array.rds"))
truth_male_3d = readRDS(file = file.path(dir.r, "male_prefecture_rate_array.rds"))

# data_set_array: n_age x n_state x n_year
# fh: a forecast horizon
# fmethod: forecasting method, ets or auto.arima

FM_TWA_FFM_fun <- function(data_set_array, fh, fmethod, Fixed_TWA)
{
  n_age = dim(data_set_array)[1]
  n_state = dim(data_set_array)[2]
  n_year = dim(data_set_array)[3]
  
  data_set_mat = NULL
  for(ik in 1:n_state)
  {
    data_set_mat = cbind(data_set_mat, data_set_array[,ik,])
    rm(ik)
  }
  
  TWA_FM_smooth_resid <- data_set_mat
  TWA_FM_smooth_resid_array = array(NA, dim = c(n_age, n_year, n_state), dimnames = list(1:n_age, 1:n_year, 1:n_state))
  for(iw in 1:n_state)
  {
    TWA_FM_smooth_resid_array[,,iw] = (TWA_FM_smooth_resid)[,((iw-1)*n_year+1):(iw*n_year)]
    rm(iw)
  }
  
  HDFTS_factor_decomp_female <- HDFTS_factor_decomp(data = aperm(TWA_FM_smooth_resid_array, c(2, 3, 1)))
  HDFTS_factor_decomp_female_factor_loading <- HDFTS_factor_decomp_female$factor_loading
  HDFTS_factor_decomp_female_factors <- HDFTS_factor_decomp_female$factors
  HDFTS_factor_decomp_female_residuals <- HDFTS_factor_decomp_female$error
  
  # forecasting
  
  n_factor <- ncol(HDFTS_factor_decomp_female_factors)
  HDFTS_factor_decomp_female_factors_forecast = matrix(NA, 1, n_factor)
  for(ik in 1:n_factor)
  {
    if(fmethod == "ets")
    {
      HDFTS_factor_decomp_female_factors_forecast[,ik] = forecast(ets(HDFTS_factor_decomp_female_factors[,ik]), h = fh)$mean[fh]
    }
    else if(fmethod == "arima")
    {
      HDFTS_factor_decomp_female_factors_forecast[,ik] = forecast(auto.arima(HDFTS_factor_decomp_female_factors[,ik]), h = fh)$mean[fh]
    }
    rm(ik)
  }
  
  HDFTS_factor_decomp_female_forecast <- matrix(NA, n_age, n_state)
  for(ik in 1:n_state)
  {
    HDFTS_factor_decomp_female_forecast[,ik] = HDFTS_factor_decomp_female_factors_forecast %*% HDFTS_factor_decomp_female_factor_loading[,,ik]
    rm(ik)
  }
  FM_FFM_forecast = HDFTS_factor_decomp_female_forecast + Fixed_TWA
  colnames(FM_FFM_forecast) = 1:n_state
  rownames(FM_FFM_forecast) = 0:(n_age - 1)
  return(FM_FFM_forecast)
}



# raw_data: mortality rate in the original scale
# smooth_data: smoothed mortality rate
# fh: forecast horizon
# fore_method: forecasting method, ets or auto.arima
# alpha: level of significance


FM_TWA_FFM_fun_int_eval_sequential <- function(raw_data, smooth_data, fh, fore_method, alpha,Fixed_TWA)
{
  n_age = dim(raw_data)[1]
  n_state = dim(raw_data)[2]
  n_year = dim(raw_data)[3]
  
  # define the test set
  
  n_validation = n_year - 11
  
  int_score_array = array(NA, dim = c((11 - fh), 3, n_state), dimnames = list(1:(11 - fh), c("ECP", "CPD", "score"), 1:n_state))
  q_alpha_mat = array(NA, dim = c(n_age, (11 - fh), n_state))
  for(iwk in 1:(11 - fh))
  {
    validation_set = 2:(n_validation + iwk - fh)
    fore_validation = array(NA, dim = c(n_age, n_state, length(validation_set)))
    for(ij in validation_set)
    {
      fore_validation[,,(ij - 1)] = FM_TWA_FFM_fun(data_set_array = (smooth_data[,,1:ij]), fh = fh, fmethod = fore_method,Fixed_TWA)
      rm(ij)
    }
    
    # compute absolute of residuals
    
    resi = abs(raw_data[,,(2 + fh):((n_validation + iwk))] - exp(fore_validation))
    q_alpha = matrix(NA, n_age, n_state)
    for(iw in 1:n_state)
    {
      for(ij in 1:n_age)
      {
        # (1 - alpha) quantile of the absolute residuals
        
        rq_fit <- rq(resi[ij,iw,] ~ 1, tau = (1 - alpha))
        
        # use information criterion to select AR order
        
        AR_p <- max(ar(resi[ij,iw,])$order, 1)
        lastlags <- tail(x = resi[ij,iw,], AR_p)
        newdf <- as.data.frame(t(lastlags))
        names(newdf) <- paste0("lag", 1:AR_p)
        
        q_alpha[ij,iw] = predict(rq_fit, newdf)
        rm(ij)
      }
      rm(iw)
    }
    colnames(q_alpha) = 1:n_state
    rownames(q_alpha) = 1:n_age
    
    point_forecast <- exp(FM_TWA_FFM_fun(data_set_array = (smooth_data[,,1:(n_validation + iwk)]), fh = fh, fmethod = fore_method,Fixed_TWA))
    for(iw in 1:n_state)
    {
      int_score_array[iwk,,iw] <- interval_score(holdout = raw_data[,iw,(n_validation + iwk + fh)], 
                                                 lb = (point_forecast - q_alpha)[,iw],
                                                 ub = (point_forecast + q_alpha)[,iw], 
                                                 alpha = alpha)
    }
    q_alpha_mat[,iwk,] = q_alpha
    rm(iwk)
  }
  return(list(q_alpha_array = q_alpha_mat, int_score_array = int_score_array))
}


### alpha = 0.2

## arima

# female

sequential_int_score_list_female_arima = sequential_int_score_list_male_arima = list()
for(iw in 1:10)
{
  sequential_int_score_list_female_arima[[iw]] = FM_TWA_FFM_fun_int_eval_sequential(raw_data = truth_female_3d, 
                                                                                smooth_data = TWA_res_female, 
                                                                                fh = iw, fore_method = "arima", alpha = 0.2,Fixed_TWA=TWA_Fixed_female)
  
  sequential_int_score_list_male_arima[[iw]] = FM_TWA_FFM_fun_int_eval_sequential(raw_data = truth_male_3d, 
                                                                              smooth_data = TWA_res_male, 
                                                                              fh = iw, fore_method = "arima", alpha = 0.2,Fixed_TWA=TWA_Fixed_male)
  print(iw); rm(iw)
}

sequential_int_score_list_female_arima_mat = sequential_int_score_list_male_arima_mat = matrix(NA, 10, 3)
for(iw in 1:10)
{
  sequential_int_score_list_female_arima_mat[iw,] = apply(sequential_int_score_list_female_arima[[iw]]$int_score_array, 2, mean)
  sequential_int_score_list_male_arima_mat[iw,]   = apply(sequential_int_score_list_male_arima[[iw]]$int_score_array, 2, mean)
  rm(iw)
}
rownames(sequential_int_score_list_female_arima_mat) = rownames(sequential_int_score_list_male_arima_mat) = 1:10
colnames(sequential_int_score_list_female_arima_mat) = colnames(sequential_int_score_list_male_arima_mat) = c("ECP", "CPD", "score")

## ets

# female

sequential_int_score_list_female_ets = sequential_int_score_list_male_ets = list()
for(iw in 1:10)
{
  sequential_int_score_list_female_ets[[iw]] = FM_TWA_FFM_fun_int_eval_sequential(raw_data = truth_female_3d, 
                                                                              smooth_data = TWA_res_female, 
                                                                              fh = iw, fore_method = "ets", alpha = 0.2,Fixed_TWA=TWA_Fixed_female)
  
  sequential_int_score_list_male_ets[[iw]] = FM_TWA_FFM_fun_int_eval_sequential(raw_data = truth_male_3d, 
                                                                            smooth_data = TWA_res_male, 
                                                                            fh = iw, fore_method = "ets", alpha = 0.2,Fixed_TWA=TWA_Fixed_male)
  print(iw); rm(iw)
}

sequential_int_score_list_female_ets_mat = sequential_int_score_list_male_ets_mat = matrix(NA, 10, 3)
for(iw in 1:10)
{
  sequential_int_score_list_female_ets_mat[iw,] = apply(sequential_int_score_list_female_ets[[iw]]$int_score_array, 2, mean)
  sequential_int_score_list_male_ets_mat[iw,]   = apply(sequential_int_score_list_male_ets[[iw]]$int_score_array, 2, mean)
  rm(iw)
}
rownames(sequential_int_score_list_female_ets_mat) = rownames(sequential_int_score_list_male_ets_mat) = 1:10
colnames(sequential_int_score_list_female_ets_mat) = colnames(sequential_int_score_list_male_ets_mat) = c("ECP", "CPD", "score")


### alpha = 0.05

## arima

# female

sequential_int_score_list_female_arima_PI_95 = sequential_int_score_list_male_arima_PI_95 = list()
for(iw in 1:10)
{
  sequential_int_score_list_female_arima_PI_95[[iw]] = FM_TWA_FFM_fun_int_eval_sequential(raw_data = truth_female_3d, 
                                                                                      smooth_data = TWA_res_female, 
                                                                                      fh = iw, fore_method = "arima", alpha = 0.05,Fixed_TWA=TWA_Fixed_female)
  
  sequential_int_score_list_male_arima_PI_95[[iw]] = FM_TWA_FFM_fun_int_eval_sequential(raw_data = truth_male_3d, 
                                                                                    smooth_data = TWA_res_male, 
                                                                                    fh = iw, fore_method = "arima", alpha = 0.05,Fixed_TWA=TWA_Fixed_male)
  print(iw); rm(iw)
}

sequential_int_score_list_female_arima_PI_95_value = sequential_int_score_list_male_arima_PI_95_value = array(NA, dim = c(47, 3, 10))
for(iw in 1:10)
{
  sequential_int_score_list_female_arima_PI_95_value[,,iw] = t(apply(sequential_int_score_list_female_arima_PI_95[[iw]]$int_score_array, c(2, 3), mean))
  sequential_int_score_list_male_arima_PI_95_value[,,iw]   = t(apply(sequential_int_score_list_male_arima_PI_95[[iw]]$int_score_array, c(2, 3), mean))
  rm(iw)
}

sequential_int_score_list_female_arima_PI_95_mat = sequential_int_score_list_male_arima_PI_95_mat = matrix(NA, 10, 3)
for(iw in 1:10)
{
  sequential_int_score_list_female_arima_PI_95_mat[iw,] = apply(sequential_int_score_list_female_arima_PI_95[[iw]]$int_score_array, 2, mean)
  sequential_int_score_list_male_arima_PI_95_mat[iw,]   = apply(sequential_int_score_list_male_arima_PI_95[[iw]]$int_score_array, 2, mean)
  rm(iw)
}
rownames(sequential_int_score_list_female_arima_PI_95_mat) = rownames(sequential_int_score_list_male_arima_PI_95_mat) = 1:10
colnames(sequential_int_score_list_female_arima_PI_95_mat) = colnames(sequential_int_score_list_male_arima_PI_95_mat) = c("ECP", "CPD", "score")

## ets

# female

sequential_int_score_list_female_ets_PI_95 = sequential_int_score_list_male_ets_PI_95 = list()
for(iw in 1:10)
{
  sequential_int_score_list_female_ets_PI_95[[iw]] = FM_TWA_FFM_fun_int_eval_sequential(raw_data = truth_female_3d, 
                                                                                    smooth_data = TWA_res_female, 
                                                                                    fh = iw, fore_method = "ets", alpha = 0.05,Fixed_TWA=TWA_Fixed_female)
  
  sequential_int_score_list_male_ets_PI_95[[iw]] = FM_TWA_FFM_fun_int_eval_sequential(raw_data = truth_male_3d, 
                                                                                  smooth_data = TWA_res_male, 
                                                                                  fh = iw, fore_method = "ets", alpha = 0.05,Fixed_TWA=TWA_Fixed_male)
  print(iw); rm(iw)
}

sequential_int_score_list_female_ets_PI_95_value = sequential_int_score_list_male_ets_PI_95_value = array(NA, dim = c(47, 3, 10))
for(iw in 1:10)
{
  sequential_int_score_list_female_ets_PI_95_value[,,iw] = t(apply(sequential_int_score_list_female_ets_PI_95[[iw]]$int_score_array, c(2, 3), mean))
  sequential_int_score_list_male_ets_PI_95_value[,,iw]   = t(apply(sequential_int_score_list_male_ets_PI_95[[iw]]$int_score_array, c(2, 3), mean))
  rm(iw)
}

sequential_int_score_list_female_ets_PI_95_mat = sequential_int_score_list_male_ets_PI_95_mat = matrix(NA, 10, 3)
for(iw in 1:10)
{
  sequential_int_score_list_female_ets_PI_95_mat[iw,] = apply(sequential_int_score_list_female_ets_PI_95[[iw]]$int_score_array, 2, mean)
  sequential_int_score_list_male_ets_PI_95_mat[iw,]   = apply(sequential_int_score_list_male_ets_PI_95[[iw]]$int_score_array, 2, mean)
  rm(iw)
}
rownames(sequential_int_score_list_female_ets_PI_95_mat) = rownames(sequential_int_score_list_male_ets_PI_95_mat) = 1:10
colnames(sequential_int_score_list_female_ets_PI_95_mat) = colnames(sequential_int_score_list_male_ets_PI_95_mat) = c("ECP", "CPD", "score")

