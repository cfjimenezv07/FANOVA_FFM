###########################################################################
# One way functional ANOVA
###########################################################################

###########################################################################
# Grand mean function
###########################################################################
#' Description for the function mu_hat_oneway
# 1. Computes the functional grand mean for one-way FANOVA.
#'
#' Parameters
#' @param data_pop1 A matrix of functional data with dimension n_age by
#'        (n_year * n_prefectures).
#' @param n_year Number of years observed for each prefecture.
#' @param n_prefectures Number of prefectures.
#' @param n_age Number of grid points (ages).
#'
#' Returns
#' @return A vector of length n_age representing the grand mean function.
###########################################################################

mu_hat_oneway <- function(data_pop1,
                          n_year,
                          n_prefectures,
                          n_age) {
  
  Mu_hat <- rep(0, n_age)
  
  for (i in 1:n_prefectures) {
    block_i <- t(
      data_pop1[, (n_year * i - (n_year - 1)):(n_year * i)]
    )
    for (k in 1:n_year) {
      Mu_hat <- Mu_hat + block_i[k, ]
    }
  }
  
  Mu_hat <- Mu_hat / (n_year * n_prefectures)
  return(Mu_hat)
}

###########################################################################
# Prefecture (row) effects
###########################################################################
#' Description for the function FRE_means_oneway
# 1. Computes the functional prefecture effects for one-way FANOVA.
#'
#' Parameters
#' @param data_pop1 Functional data matrix (n_age by n_year * n_prefectures).
#' @param n_year Number of years observed per prefecture.
#' @param n_prefectures Number of prefectures.
#' @param n_age Number of grid points (ages).
#'
#' Returns
#' @return A matrix of dimension n_age by n_prefectures containing
#'         the functional prefecture effects.
###########################################################################

FRE_means_oneway <- function(data_pop1,
                             n_year,
                             n_prefectures,
                             n_age) {
  
  Mu_hat <- mu_hat_oneway(
    data_pop1,
    n_year,
    n_prefectures,
    n_age
  )
  
  FRE_means <- matrix(0, nrow = n_age, ncol = n_prefectures)
  
  for (i in 1:n_prefectures) {
    
    alpha_hat <- rep(0, n_age)
    
    block_i <- t(
      data_pop1[, (n_year * i - (n_year - 1)):(n_year * i)]
    )
    
    for (k in 1:n_year) {
      alpha_hat <- alpha_hat + block_i[k, ]
    }
    
    FRE_means[, i] <- (alpha_hat / n_year) - Mu_hat
  }
  
  return(FRE_means)
}

###########################################################################
# One way FANOVA deterministic components
###########################################################################
#' Description for the function FANOVA_oneway
# 1. Computes the deterministic part of the one-way functional ANOVA.
# 2. Returns the grand mean, prefecture effects, and their sum.
#'
#' Parameters
#' @param data_pop1 Functional data matrix (n_age by n_year * n_prefectures).
#' @param year Vector of observed years.
#' @param age Vector of ages.
#' @param n_prefectures Number of prefectures.
#'
#' Returns
#' @return A list containing:
#'   \item{GE_mean}{Grand mean function}
#'   \item{FRE_mean}{Prefecture effects (n_prefectures by n_age)}
#'   \item{Deterministic}{Matrix representing mu(a) + alpha_i(a)}
###########################################################################

FANOVA_oneway <- function(data_pop1,
                          year = 1959:2020,
                          age = 0:100,
                          n_prefectures = 51) {
  
  n_year <- length(year)
  n_age  <- length(age)
  
  Mu_hat <- mu_hat_oneway(
    data_pop1,
    n_year,
    n_prefectures,
    n_age
  )
  
  FRE_mean <- FRE_means_oneway(
    data_pop1,
    n_year,
    n_prefectures,
    n_age
  )
  
  Deterministic <- matrix(0,
                          nrow = n_age,
                          ncol = n_year * n_prefectures)
  
  for (i in 1:n_prefectures) {
    cols_i <- (n_year * i - (n_year - 1)):(n_year * i)
    for (k in 1:n_year) {
      Deterministic[, cols_i[k]] <- Mu_hat + FRE_mean[, i]
    }
  }
  
  return(list(
    GE_mean = Mu_hat,
    FRE_mean = t(FRE_mean),
    Deterministic = Deterministic
  ))
}


###########################################################################
# One way functional residuals
###########################################################################
#' Description for the function One_way_Residuals
# 1. Computes the functional residuals from the one-way FANOVA.
# 2. Confirms that mu(a) + alpha_i(a) + epsilon_ik(a) reconstructs the data.
#'
#' Parameters
#' @param data_pop1 Functional data matrix (n_age by n_year * n_prefectures).
#' @param n_prefectures Number of prefectures.
#' @param n_year Number of years observed.
#' @param n_age Number of grid points (ages).
#'
#' Returns
#' @return A list containing:
#'   \item{Residuals}{Functional residuals with same dimension as data_pop1}
#'   \item{Reconstructed_Data}{mu + alpha_i + epsilon_ik}
#'   \item{Reconstruction_OK}{Logical confirmation of exact reconstruction}
###########################################################################

One_way_Residuals <- function(data_pop1,
                              n_prefectures,
                              n_year,
                              n_age) {
  
  Mu_hat <- mu_hat_oneway(
    data_pop1,
    n_year,
    n_prefectures,
    n_age
  )
  
  FRE_mean <- FRE_means_oneway(
    data_pop1,
    n_year,
    n_prefectures,
    n_age
  )
  
  Residuals <- matrix(0,
                      nrow = n_age,
                      ncol = n_year * n_prefectures)
  
  Reconstructed_Data <- matrix(0,
                               nrow = n_age,
                               ncol = n_year * n_prefectures)
  
  for (i in 1:n_prefectures) {
    
    cols_i <- (n_year * i - (n_year - 1)):(n_year * i)
    
    for (k in 1:n_year) {
      
      Residuals[, cols_i[k]] <-
        data_pop1[, cols_i[k]] - Mu_hat - FRE_mean[, i]
      
      Reconstructed_Data[, cols_i[k]] <-
        Mu_hat + FRE_mean[, i] + Residuals[, cols_i[k]]
    }
  }
  
  Reconstruction_OK <- isTRUE(
    all.equal(Reconstructed_Data, data_pop1)
  )
  
  return(list(
    Residuals = Residuals,
    Reconstructed_Data = Reconstructed_Data,
    Reconstruction_OK = Reconstruction_OK
  ))
}



