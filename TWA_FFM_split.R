

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

# truth_3d: [96 x 47 x 49] (Age x Year x Prefecture)
# female_smooth / male_smooth: [96 x 2303] (Ages x (Years*Prefectures))
female_prefecture_rate_combo_smooth   <- readRDS(file.path(dir.r, "female_prefecture_rate_combo_smooth.rds"))
male_prefecture_rate_combo_smooth     <- readRDS(file.path(dir.r, "male_prefecture_rate_combo_smooth.rds"))

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


# data_set_array: n_age x n_state x n_year
# fh: a forecast horizon
# fmethod: forecasting method, ets or auto.arima

# raw_data = truth_female_3d
# smooth_data = TWA_res_female
# iw=1
# fh = iw
# fmethod =fore_method = "arima"
# alpha = 0.2
# Fixed_TWA=TWA_Fixed_female

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

# female_prefecture_rate_array_smooth_FM_FFM <- FM_TWA_FFM_fun(data_set_array = TWA_res_female, fh = 10, fmethod = "ets",Fixed_TWA=TWA_Fixed_female)
# male_prefecture_rate_array_smooth_FM_FFM   <- FM_TWA_FFM_fun(data_set_array = TWA_res_male, fh = 10, fmethod = "ets",Fixed_TWA=TWA_Fixed_male)

##############
# calibration
##############

tune_para_select <- function(resi, tune_para, nominal, resi_sd)
{
  n_age = length(resi_sd)
  ind = matrix(NA, n_age, ncol(resi))
  for(ijk in 1:ncol(resi))
  {
    ind[,ijk] = ifelse(between(x = resi[,ijk], left = -tune_para * resi_sd, right = tune_para * resi_sd), 1, 0)
    rm(ijk)
  }
  return(abs(sum(ind)/(n_age * ncol(resi)) - nominal))
}

# training_data: 1:28
# validation_data: 29:39
# test_data: 40:49

# raw_data: mortality rate in the original scale
# smooth_data: smoothed mortality rate
# fh: forecast horizon
# fore_method: forecasting method, ets or auto.arima
# alpha: level of significance



FM_TWA_FFM_fun_int_eval <- function(raw_data, smooth_data, fh, fore_method, alpha,Fixed_TWA)
{
  n_age = dim(raw_data)[1]
  n_state = dim(raw_data)[2]
  n_year = dim(raw_data)[3]

  # validation set

  smooth_FM_FFM_mat = array(NA, dim = c(n_age, n_state, (12 - fh)),
                            dimnames = list(1:n_age, 1:n_state, 1:(12 - fh)))
  for(iw in 1:(12 - fh))
  {
    smooth_FM_FFM_mat[,,iw] <- FM_TWA_FFM_fun(data_set_array = smooth_data[,,1:((n_year - 22) + iw)], fh = fh, fmethod = fore_method,Fixed_TWA)
    rm(iw)
  }
  point_forecast <- exp(smooth_FM_FFM_mat)
  validation_set <- array(raw_data[,,(n_year - (21 - fh)):(n_year - 10)], dim = c(n_age, n_state, (12 - fh)))
  FM_FFM_err <- validation_set - point_forecast

  # compute pointwise sd

  FM_FFM_quantile = FM_FFM_sd = matrix(NA, n_age, n_state)
  for(iw in 1:n_state)
  {
    FM_FFM_quantile[,iw] = apply(abs(FM_FFM_err[,iw,]), 1, quantile, (1 - alpha))
    FM_FFM_sd[,iw] = apply(FM_FFM_err[,iw,], 1, sd)
    rm(iw)
  }
  colnames(FM_FFM_quantile) = colnames(FM_FFM_sd) = 1:n_state
  rownames(FM_FFM_quantile) = rownames(FM_FFM_sd) = 1:n_age

  # try several optimisers

  tune_para_value = vector("numeric", n_state)
  for(iw in 1:n_state)
  {
    obj_BFGS = optim(par = 1, fn = tune_para_select, method = "BFGS", resi = FM_FFM_err[,iw,],
                     nominal = 1 - alpha, resi_sd = FM_FFM_sd[,iw])
    obj_NM = optim(par = 1, fn = tune_para_select, method = "Nelder-Mead", resi = FM_FFM_err[,iw,],
                   nominal = 1 - alpha, resi_sd = FM_FFM_sd[,iw])
    obj_Brent = optim(par = 1, fn = tune_para_select, method = "Brent", resi = FM_FFM_err[,iw,],
                      nominal = 1 - alpha, lower = 0, upper = 10, resi_sd = FM_FFM_sd[,iw])
    tune_para_obj = min(c(obj_BFGS$value, obj_NM$value, obj_Brent$value))
    tune_para_value[iw] = c(obj_BFGS$par, obj_NM$par, obj_Brent$par)[which.min(c(obj_BFGS$value, obj_NM$value, obj_Brent$value))]
    rm(iw)
  }
  rm(point_forecast); rm(validation_set); rm(FM_FFM_err)

  # forecasting set

  smooth_FM_FFM_mat_fore = array(NA, dim = c(n_age, n_state, (11 - fh)),
                                 dimnames = list(1:n_age, 1:n_state, 1:(11 - fh)))
  for(iw in 1:(11 - fh))
  {
    smooth_FM_FFM_mat_fore[,,iw] <- FM_TWA_FFM_fun(data_set_array = smooth_data[,,1:((n_year - 11) + iw)], fh = fh, fmethod = fore_method,Fixed_TWA)
    rm(iw)
  }

  point_forecast <- exp(smooth_FM_FFM_mat_fore)
  lb_array_sd = ub_array_sd = lb_array_quantile = ub_array_quantile = array(NA, dim = c(n_age, (11 - fh), n_state), dimnames = list(1:n_age, 1:(11 - fh), 1:n_state))
  int_score_sd = int_score_quantile = matrix(NA, n_state, 3)
  for(iw in 1:n_state)
  {
    lb_array_sd[,,iw] = point_forecast[,iw,] - tune_para_value[iw] * FM_FFM_sd[,iw]
    ub_array_sd[,,iw] = point_forecast[,iw,] + tune_para_value[iw] * FM_FFM_sd[,iw]

    lb_array_quantile[,,iw] = point_forecast[,iw,] - FM_FFM_quantile[,iw]
    ub_array_quantile[,,iw] = point_forecast[,iw,] + FM_FFM_quantile[,iw]

    int_score_sd[iw,] = interval_score(holdout = raw_data[,iw,(n_year - (10 - fh)):n_year], lb = lb_array_sd[,,iw],
                                       ub = ub_array_sd[,,iw], alpha = alpha)
    int_score_quantile[iw,] = interval_score(holdout = raw_data[,iw,(n_year - (10 - fh)):n_year], lb = lb_array_quantile[,,iw],
                                             ub = ub_array_quantile[,,iw], alpha = alpha)
    rm(iw)
  }
  rownames(int_score_sd) = rownames(int_score_quantile) = 1:n_state
  colnames(int_score_sd) = colnames(int_score_quantile) = c("ECP", "CPD", "score")
  return(list(int_score_sd = int_score_sd, int_score_quantile = int_score_quantile))
}

################
### alpha = 0.2
################
## female

# arima

FM_TWA_FFM_fun_int_score_female_array_arima_sd = FM_TWA_FFM_fun_int_score_female_array_arima_quantile = array(NA, dim = c(47, 3, 10), dimnames = list(1:47, c("ECP", "CPD", "score"), 1:10))
for(iw in 1:10)
{
  dum = FM_TWA_FFM_fun_int_eval(raw_data = truth_female_3d,
                            smooth_data = TWA_res_female,
                            fh = iw, fore_method = "arima", alpha = 0.2,Fixed_TWA=TWA_Fixed_female)
  FM_TWA_FFM_fun_int_score_female_array_arima_sd[,,iw] = dum$int_score_sd
  FM_TWA_FFM_fun_int_score_female_array_arima_quantile[,,iw] = dum$int_score_quantile
  print(iw); rm(iw); rm(dum)
}

FM_TWA_FFM_fun_int_score_female_array_arima_sd_mean = t(apply(FM_TWA_FFM_fun_int_score_female_array_arima_sd, c(2, 3), mean))
FM_TWA_FFM_fun_int_score_female_array_arima_quantile_mean = t(apply(FM_TWA_FFM_fun_int_score_female_array_arima_quantile, c(2, 3), mean))

# ets

FM_TWA_FFM_fun_int_score_female_array_ets_sd = FM_TWA_FFM_fun_int_score_female_array_ets_quantile = array(NA, dim = c(47, 3, 10), dimnames = list(1:47, c("ECP", "CPD", "score"), 1:10))
for(iw in 1:10)
{
  dum = FM_TWA_FFM_fun_int_eval(raw_data = truth_female_3d,
                            smooth_data = TWA_res_female,
                            fh = iw, fore_method = "ets", alpha = 0.2,Fixed_TWA=TWA_Fixed_female)
  FM_TWA_FFM_fun_int_score_female_array_ets_sd[,,iw] = dum$int_score_sd
  FM_TWA_FFM_fun_int_score_female_array_ets_quantile[,,iw] = dum$int_score_quantile
  print(iw); rm(iw); rm(dum)
}

FM_TWA_FFM_fun_int_score_female_array_ets_sd_mean = t(apply(FM_TWA_FFM_fun_int_score_female_array_ets_sd, c(2, 3), mean))
FM_TWA_FFM_fun_int_score_female_array_ets_quantile_mean = t(apply(FM_TWA_FFM_fun_int_score_female_array_ets_quantile, c(2, 3), mean))

## male

# arima

FM_TWA_FFM_fun_int_score_male_array_arima_sd = FM_TWA_FFM_fun_int_score_male_array_arima_quantile = array(NA, dim = c(47, 3, 10), dimnames = list(1:47, c("ECP", "CPD", "score"), 1:10))
for(iw in 1:10)
{
  dum = FM_TWA_FFM_fun_int_eval(raw_data = truth_male_3d,
                            smooth_data = TWA_res_male,
                            fh = iw, fore_method = "arima", alpha = 0.2,Fixed_TWA=TWA_Fixed_male)
  FM_TWA_FFM_fun_int_score_male_array_arima_sd[,,iw] = dum$int_score_sd
  FM_TWA_FFM_fun_int_score_male_array_arima_quantile[,,iw] = dum$int_score_quantile
  print(iw); rm(iw); rm(dum)
}

FM_TWA_FFM_fun_int_score_male_array_arima_sd_mean = t(apply(FM_TWA_FFM_fun_int_score_male_array_arima_sd, c(2, 3), mean))
FM_TWA_FFM_fun_int_score_male_array_arima_quantile_mean = t(apply(FM_TWA_FFM_fun_int_score_male_array_arima_quantile, c(2, 3), mean))

# ets

FM_TWA_FFM_fun_int_score_male_array_ets_sd = FM_TWA_FFM_fun_int_score_male_array_ets_quantile = array(NA, dim = c(47, 3, 10), dimnames = list(1:47, c("ECP", "CPD", "score"), 1:10))
for(iw in 1:10)
{
  dum = FM_TWA_FFM_fun_int_eval(raw_data = truth_male_3d,
                            smooth_data = TWA_res_male,
                            fh = iw, fore_method = "ets", alpha = 0.2,Fixed_TWA=TWA_Fixed_male)
  FM_TWA_FFM_fun_int_score_male_array_ets_sd[,,iw] = dum$int_score_sd
  FM_TWA_FFM_fun_int_score_male_array_ets_quantile[,,iw] = dum$int_score_quantile
  print(iw); rm(iw); rm(dum)
}

FM_TWA_FFM_fun_int_score_male_array_ets_sd_mean = t(apply(FM_TWA_FFM_fun_int_score_male_array_ets_sd, c(2, 3), mean))
FM_TWA_FFM_fun_int_score_male_array_ets_quantile_mean = t(apply(FM_TWA_FFM_fun_int_score_male_array_ets_quantile, c(2, 3), mean))

#################
### alpha = 0.05
#################

## female

# arima

FM_TWA_FFM_fun_int_score_female_array_arima_sd_alpha_0.05 = FM_TWA_FFM_fun_int_score_female_array_arima_quantile_alpha_0.05 = array(NA, dim = c(47, 3, 10), dimnames = list(1:47, c("ECP", "CPD", "score"), 1:10))
for(iw in 1:10)
{
  dum = FM_TWA_FFM_fun_int_eval(raw_data = truth_female_3d,
                            smooth_data = TWA_res_female,
                            fh = iw, fore_method = "arima", alpha = 0.05,Fixed_TWA=TWA_Fixed_female)
  FM_TWA_FFM_fun_int_score_female_array_arima_sd_alpha_0.05[,,iw] = dum$int_score_sd
  FM_TWA_FFM_fun_int_score_female_array_arima_quantile_alpha_0.05[,,iw] = dum$int_score_quantile
  print(iw); rm(iw); rm(dum)
}

FM_TWA_FFM_fun_int_score_female_array_arima_sd_alpha_0.05_mean = t(apply(FM_TWA_FFM_fun_int_score_female_array_arima_sd_alpha_0.05, c(2, 3), mean))
FM_TWA_FFM_fun_int_score_female_array_arima_quantile_alpha_0.05_mean = t(apply(FM_TWA_FFM_fun_int_score_female_array_arima_quantile_alpha_0.05, c(2, 3), mean))

# ets

FM_TWA_FFM_fun_int_score_female_array_ets_sd_alpha_0.05 = FM_TWA_FFM_fun_int_score_female_array_ets_quantile_alpha_0.05 = array(NA, dim = c(47, 3, 10), dimnames = list(1:47, c("ECP", "CPD", "score"), 1:10))
for(iw in 1:10)
{
  dum = FM_TWA_FFM_fun_int_eval(raw_data = truth_female_3d,
                            smooth_data = TWA_res_female,
                            fh = iw, fore_method = "ets", alpha = 0.05,Fixed_TWA=TWA_Fixed_female)
  FM_TWA_FFM_fun_int_score_female_array_ets_sd_alpha_0.05[,,iw] = dum$int_score_sd
  FM_TWA_FFM_fun_int_score_female_array_ets_quantile_alpha_0.05[,,iw] = dum$int_score_quantile
  print(iw); rm(iw); rm(dum)
}

FM_TWA_FFM_fun_int_score_female_array_ets_sd_alpha_0.05_mean = t(apply(FM_TWA_FFM_fun_int_score_female_array_ets_sd_alpha_0.05, c(2, 3), mean))
FM_TWA_FFM_fun_int_score_female_array_ets_quantile_alpha_0.05_mean = t(apply(FM_TWA_FFM_fun_int_score_female_array_ets_quantile_alpha_0.05, c(2, 3), mean))

## male

# arima

FM_TWA_FFM_fun_int_score_male_array_arima_sd_alpha_0.05 = FM_TWA_FFM_fun_int_score_male_array_arima_quantile_alpha_0.05 = array(NA, dim = c(47, 3, 10), dimnames = list(1:47, c("ECP", "CPD", "score"), 1:10))
for(iw in 1:10)
{
  dum = FM_TWA_FFM_fun_int_eval(raw_data = truth_male_3d,
                            smooth_data = TWA_res_male,
                            fh = iw, fore_method = "arima", alpha = 0.05,Fixed_TWA=TWA_Fixed_male)
  FM_TWA_FFM_fun_int_score_male_array_arima_sd_alpha_0.05[,,iw] = dum$int_score_sd
  FM_TWA_FFM_fun_int_score_male_array_arima_quantile_alpha_0.05[,,iw] = dum$int_score_quantile
  print(iw); rm(iw); rm(dum)
}

FM_TWA_FFM_fun_int_score_male_array_arima_sd_alpha_0.05_mean = t(apply(FM_TWA_FFM_fun_int_score_male_array_arima_sd_alpha_0.05, c(2, 3), mean))
FM_TWA_FFM_fun_int_score_male_array_arima_quantile_alpha_0.05_mean = t(apply(FM_TWA_FFM_fun_int_score_male_array_arima_quantile_alpha_0.05, c(2, 3), mean))

# ets

FM_TWA_FFM_fun_int_score_male_array_ets_sd_alpha_0.05 = FM_TWA_FFM_fun_int_score_male_array_ets_quantile_alpha_0.05 = array(NA, dim = c(47, 3, 10), dimnames = list(1:47, c("ECP", "CPD", "score"), 1:10))
for(iw in 1:10)
{
  dum = FM_TWA_FFM_fun_int_eval(raw_data = truth_male_3d,
                            smooth_data = TWA_res_male,
                            fh = iw, fore_method = "ets", alpha = 0.05,Fixed_TWA=TWA_Fixed_male)
  FM_TWA_FFM_fun_int_score_male_array_ets_sd_alpha_0.05[,,iw] = dum$int_score_sd
  FM_TWA_FFM_fun_int_score_male_array_ets_quantile_alpha_0.05[,,iw] = dum$int_score_quantile
  print(iw); rm(iw); rm(dum)
}

FM_TWA_FFM_fun_int_score_male_array_ets_sd_alpha_0.05_mean = t(apply(FM_TWA_FFM_fun_int_score_male_array_ets_sd_alpha_0.05, c(2, 3), mean))
FM_TWA_FFM_fun_int_score_male_array_ets_quantile_alpha_0.05_mean = t(apply(FM_TWA_FFM_fun_int_score_male_array_ets_quantile_alpha_0.05, c(2, 3), mean))

