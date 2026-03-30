# TWA+ FTSM CFJV model (Custom ARIMA Bootstrap Version)

################################################################################
# 0. CUSTOM FORECAST FUNCTIONS & NAMESPACE SETUP
################################################################################
library(forecast)

# Accessing internal forecast functions
getxreg         <- forecast:::getxreg
arima.string    <- forecast:::arima.string
fitted.Arima    <- forecast:::fitted.Arima
future_msts     <- forecast:::future_msts
copy_msts       <- forecast:::copy_msts
residuals.Arima <- forecast:::residuals.Arima

# Your custom forecast_Arima implementation
forecast_Arima <- function (object, h = ifelse(object$arma[5] > 1, 2 * object$arma[5], 10), 
                            level = c(80, 95), fan = FALSE, xreg = NULL, lambda = object$lambda, 
                            bootstrap = FALSE, npaths = 5000, biasadj = NULL, ...) 
{
  all.args <- names(formals())
  user.args <- names(match.call())[-1L]
  check <- user.args %in% all.args
  if (!all(check)) {
    error.args <- user.args[!check]
    warning(sprintf("The non-existent %s arguments will be ignored.", error.args))
  }
  use.drift <- is.element("drift", names(object$coef))
  x <- object$x <- forecast::getResponse(object)
  usexreg <- (use.drift | is.element("xreg", names(object)))
  if (!is.null(xreg) && usexreg) {
    if (!is.numeric(xreg)) stop("xreg should be a numeric matrix or a numeric vector")
    xreg <- as.matrix(xreg)
    if (is.null(colnames(xreg))) {
      colnames(xreg) <- if (ncol(xreg) == 1) "xreg" else paste("xreg", 1:ncol(xreg), sep = "")
    }
    origxreg <- xreg <- as.matrix(xreg)
    h <- nrow(xreg)
  } else {
    if (!is.null(xreg)) {
      warning("xreg not required by this model, ignoring the provided regressors")
      xreg <- NULL
    }
    origxreg <- NULL
  }
  if (fan) {
    level <- seq(51, 99, by = 3)
  } else {
    if (min(level) > 0 & max(level) < 1) {
      level <- 100 * level
    } else if (min(level) < 0 | max(level) > 99.99) {
      stop("Confidence limit out of range")
    }
  }
  level <- sort(level)
  if (use.drift) {
    n <- length(x)
    if (!is.null(xreg)) {
      xreg <- `colnames<-`(cbind(drift = (1:h) + n, xreg), 
                           make.unique(c("drift", if (is.null(colnames(xreg)) && !is.null(xreg)) rep("", NCOL(xreg)) else colnames(xreg))))
    } else {
      xreg <- `colnames<-`(as.matrix((1:h) + n), "drift")
    }
  }
  if (!is.null(object$constant)) {
    if (object$constant) pred <- list(pred = rep(x[1], h), se = rep(0, h))
    else stop("Strange value of object$constant")
  } else if (usexreg) {
    if (is.null(xreg)) stop("No regressors provided")
    object$call$xreg <- getxreg(object)
    if (NCOL(xreg) != NCOL(object$call$xreg)) stop("Number of regressors does not match fitted model")
    if (!identical(colnames(xreg), colnames(object$call$xreg))) {
      warning("xreg contains different column names from the xreg used in training.")
    }
    pred <- predict(object, n.ahead = h, newxreg = xreg)
  } else {
    pred <- predict(object, n.ahead = h)
  }
  
  if (!is.null(x)) {
    tspx <- tsp(x)
    nx <- max(which(!is.na(x)))
    if (nx != length(x) | is.null(tsp(pred$pred)) | is.null(tsp(pred$se))) {
      tspx[2] <- time(x)[nx]
      start.f <- tspx[2] + 1/tspx[3]
      pred$pred <- ts(pred$pred, frequency = tspx[3], start = start.f)
      pred$se <- ts(pred$se, frequency = tspx[3], start = start.f)
    }
  }
  
  nint <- length(level)
  if (bootstrap) {
    sim <- matrix(NA, nrow = npaths, ncol = h)
    for (i in 1:npaths) sim[i, ] <- simulate(object, nsim = h, bootstrap = TRUE, xreg = origxreg, lambda = lambda)
    lower <- apply(sim, 2, quantile, 0.5 - level/200, type = 8)
    upper <- apply(sim, 2, quantile, 0.5 + level/200, type = 8)
    if (nint > 1L) { lower <- t(lower); upper <- t(upper) }
    else { lower <- matrix(lower, ncol = 1); upper <- matrix(upper, ncol = 1) }
  } else {
    lower <- matrix(NA, ncol = nint, nrow = length(pred$pred))
    upper <- lower
    for (i in 1:nint) {
      qq <- qnorm(0.5 * (1 + level[i]/100))
      lower[, i] <- pred$pred - qq * pred$se
      upper[, i] <- pred$pred + qq * pred$se
    }
  }
  
  colnames(lower) <- colnames(upper) <- paste(level, "%", sep = "")
  lower <- ts(lower); upper <- ts(upper)
  tsp(lower) <- tsp(upper) <- tsp(pred$pred)
  method <- arima.string(object, padding = FALSE)
  seriesname <- if (!is.null(object$series)) object$series else if (!is.null(object$call$x)) object$call$x else object$call$y
  fits <- fitted.Arima(object)
  
  if (!is.null(lambda) & is.null(object$constant)) {
    pred$pred <- InvBoxCox(pred$pred, lambda, biasadj, pred$se^2)
    if (!bootstrap) { lower <- InvBoxCox(lower, lambda); upper <- InvBoxCox(upper, lambda) }
  }
  
  return(structure(list(method = method, model = object, level = level, 
                        mean = future_msts(x, pred$pred), sim = if(bootstrap) sim else NULL, 
                        lower = future_msts(x, lower), upper = future_msts(x, upper), 
                        x = x, series = seriesname, fitted = copy_msts(x, fits), 
                        residuals = copy_msts(x, residuals.Arima(object))), class = "forecast"))
}

################################################################################
# 1. SETUP & DATA PREPARATION
################################################################################
setwd("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Rcodes")
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/FANOVA_FFM/Results"

source("1_load_packages.R")
source("FFM.R") 

state_names <- c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
                 "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", 
                 "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", 
                 "Shizuoka", "Aichi", "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", 
                 "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi", 
                 "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", 
                 "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

ages <- 0:95; years <- 1975:2023; n_years <- length(years)

# This datasets can be obtain from 2_read_data.R
female_prefecture_rate_combo_smooth <- readRDS(file.path(dir.r, "female_prefecture_rate_combo_smooth.rds"))
male_prefecture_rate_combo_smooth   <- readRDS(file.path(dir.r, "male_prefecture_rate_combo_smooth.rds"))
truth_female_3d = readRDS(file = file.path(dir.r, "female_prefecture_rate_array.rds"))
truth_male_3d = readRDS(file = file.path(dir.r, "male_prefecture_rate_array.rds"))

################################################################################
# 2. TWO-WAY ANOVA DECOMPOSITION
################################################################################
FANOVA_means <- FANOVA(log(male_prefecture_rate_combo_smooth), log(female_prefecture_rate_combo_smooth), 
                       year = years, age = ages, n_prefectures = 47, n_populations = 2)

FANOVA_residuals_raw <- Two_way_Residuals_means(log(male_prefecture_rate_combo_smooth), log(female_prefecture_rate_combo_smooth), 
                                                year = years, age = ages, n_prefectures = 47, n_populations = 2)

FANOVA_res_female_3D <- array(FANOVA_residuals_raw$residuals2_mean, dim = c(n_years, 47, 96))
FANOVA_res_male_3D   <- array(FANOVA_residuals_raw$residuals1_mean, dim = c(n_years, 47, 96))

FANOVA_deterministic_male_3D   <- array(NA, dim = c(49, 47, 96))
FANOVA_deterministic_female_3D <- array(NA, dim = c(49, 47, 96))

for(i in 1:47) {
  det_val_m <- FANOVA_means$FGE_mean + FANOVA_means$FRE_mean[i,] + FANOVA_means$FCE_mean[1,]
  det_val_f <- FANOVA_means$FGE_mean + FANOVA_means$FRE_mean[i,] + FANOVA_means$FCE_mean[2,]
  for(t in 1:49) {
    FANOVA_deterministic_male_3D[t, i, ]   <- det_val_m
    FANOVA_deterministic_female_3D[t, i, ] <- det_val_f
  }
}

################################################################################
# 3. FORECASTING ENGINE
################################################################################
library(doMC)
registerDoMC(detectCores() - 2)

select_k_evr <- function(tau = 10^-2, eigenvalue) {
  k_max = length(eigenvalue); k_all = rep(0, k_max-1)
  for(k in 1:(k_max-1)) {
    k_all[k] = (eigenvalue[k+1]/eigenvalue[k]) * ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + 
      ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
  }
  return(which.min(k_all))
}

pref_fpca_forecast <- function(residuals_3d_slice, fh, B=5000, select_K="EVR", K=6) {
  med_polish_resi = t(residuals_3d_slice)
  med_polish_resi_cov = cov(t(med_polish_resi))
  med_polish_resi_eigen = eigen(med_polish_resi_cov)
  
  retain_component = if(select_K == "Fixed") K else select_k_evr(tau = 10^-2, eigenvalue = med_polish_resi_eigen$values)
  
  basis = as.matrix(med_polish_resi_eigen$vectors[, 1:retain_component])
  score = crossprod(med_polish_resi, basis)
  
  score_fore = matrix(NA, retain_component, fh)
  for(ik in 1:retain_component) {
    # Using your requested custom function call
    dum = forecast_Arima(object = auto.arima(score[, ik]), h = fh, bootstrap = TRUE, npaths = B)
    score_fore[ik, ] = as.vector(dum$mean)
  }
  
  return(t(basis %*% score_fore)) 
}

################################################################################
# 4. EVALUATION LOOP
################################################################################
max_h <- 10; n_training_ini <- 39 
Boot_male_pref <- vector("list", 47); Boot_female_pref <- vector("list", 47)

for (i in 1:47) {
  cat(sprintf("Processing %s...\n", state_names[i]))
  res_m_slice <- FANOVA_res_male_3D[, i, ]; res_f_slice <- FANOVA_res_female_3D[, i, ]
  
  out <- foreach(iwk = 1:max_h) %dopar% {
    train_end <- n_training_ini + iwk - 1; fh_current <- n_years - train_end
    
    fore_res_m <- pref_fpca_forecast(res_m_slice[1:train_end, ], fh = fh_current)
    fore_res_f <- pref_fpca_forecast(res_f_slice[1:train_end, ], fh = fh_current)
    
    pm <- pf <- matrix(NA, fh_current, length(ages))
    for(h in 1:fh_current) {
      ty <- train_end + h
      pm[h, ] <- exp(FANOVA_deterministic_male_3D[ty, i, ] + fore_res_m[h, ])
      pf[h, ] <- exp(FANOVA_deterministic_female_3D[ty, i, ] + fore_res_f[h, ])
    }
    list(m = pm, f = pf)
  }
  Boot_male_pref[[i]] <- lapply(out, function(x) x$m)
  Boot_female_pref[[i]] <- lapply(out, function(x) x$f)
}

################################################################################
# 5. ERROR CALCULATION
################################################################################
err_m_mae <- err_f_mae <- err_m_rmse <- err_f_rmse <- lapply(1:max_h, function(x) numeric())

for (i in 1:47) {
  for (iwk in 1:max_h) {
    train_end <- n_training_ini + iwk - 1; fh_current <- n_years - train_end
    for (h in 1:fh_current) {
      target_yr <- train_end + h
      act_m <- truth_male_3d[, i, target_yr]; act_f <- truth_female_3d[, i, target_yr]
      pr_m <- Boot_male_pref[[i]][[iwk]][h, ]; pr_f <- Boot_female_pref[[i]][[iwk]][h, ]
      
      err_m_mae[[h]]  <- c(err_m_mae[[h]], mean(abs(pr_m - act_m)))
      err_f_mae[[h]]  <- c(err_f_mae[[h]], mean(abs(pr_f - act_f)))
      err_m_rmse[[h]] <- c(err_m_rmse[[h]], sqrt(mean((pr_m - act_m)^2)))
      err_f_rmse[[h]] <- c(err_f_rmse[[h]], sqrt(mean((pr_f - act_f)^2)))
    }
  }
}

final_summary <- data.frame(
  fh = as.character(1:max_h),
  Male_MAE = sapply(err_m_mae, mean), Male_RMSE = sapply(err_m_rmse, mean),
  Female_MAE = sapply(err_f_mae, mean), Female_RMSE = sapply(err_f_rmse, mean)
)

final_summary <- rbind(final_summary, 
                       data.frame(fh="Mean", t(colMeans(final_summary[,-1]))),
                       data.frame(fh="Median", t(apply(final_summary[,-1], 2, median))))

print(final_summary)
# saveRDS(final_summary, file.path(dir.r, "TWA_FFM_final_summary_arima.rds"))
