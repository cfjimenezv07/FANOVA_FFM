
##########################################################################
# novel way for selecting the number of retained factor loadings in HDFTS
##########################################################################

select_K_new <- function(eigenvalue, index, sample_size, no_pop)
{
  return((eigenvalue[index])/sample_size + index * (max(sample_size, no_pop)^(-0.5)))
}

###############################
# Numbers of factors selection
###############################

### The method of den Reijer et al. 2021 ###

k_criterion <- function(eigenvalue, plot = TRUE)
{
  n = length(eigenvalue)
  Hn = sum(1/(1:n))
  
  eigen_diff = eigenvalue[-n] - eigenvalue[-1]
  lambda_bar = 1/((2:n)*Hn)
  
  n_cross = which(eigen_diff <= lambda_bar)[1]
  
  if(n_cross > 1)
  {
    n_return = n_cross - 1
  } 
  else
  {
    n_return = 1
  }
  
  if(plot)
  {
    plot(eigen_diff[1:(n_cross+6)], type = "b", xlab = expression("Eigenvalue number" ~ italic(k)), ylab = expression("Eigenvalue" ~~~ lambda))
    lines(lambda_bar[1:(n_cross+6)], type ="b", col = 2)
    abline(v = n_return, lty = 2, col = 4)
    legend("topright", lty = c(1,1), lwd = c(2,2), col = c(1,2), c(expression(lambda ~~ "(eigenvalue)  "), expression(bar(lambda) ~~ "(threshold)  ")))
  }
  return(n_return)
}

### Eigenratio $k$ selection method ###
# tau is the threshold of eigenvalue to cut, such as 0.001

select_K <- function(tau, eigenvalue)
{
  k_max = length(eigenvalue)
  k_all = rep(0, k_max-1)
  for(k in 1:(k_max-1))
  {
    k_all[k] = (eigenvalue[k+1]/eigenvalue[k])*ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
  }
  K_hat = which.min(k_all)
  return(K_hat)
}

# novel way for selecting the number of retained factor loadings in HDFTS

select_K_new <- function(eigenvalue, index, sample_size, no_pop)
{
    return((eigenvalue[index])/sample_size + index * (max(sample_size, no_pop)^(-0.5)))
}

# Three-ways of selecting the number of factors
# data: (sample size x no_pop x no_grid)
# data: (sample size x no_pop x no_grid)

HDFTS_factor_decomp <- function(data)
{
  sample_size = dim(data)[1]
  no_pop = dim(data)[2]
  D_val = dim(data)[3]
  
  Delta = matrix(0, sample_size, sample_size)
  for(ij in 1:no_pop)
  {
    temp = matrix(NA, sample_size, sample_size)
    for(t in 1:sample_size)
    {
      for(s in 1:sample_size)
      {
        temp[t,s] = matrix(data[t,ij,], nrow = 1) %*% matrix(data[s,ij,], ncol = 1)
      }
    }
    Delta = Delta + temp/D_val
    rm(temp)
  }
  rm(t); rm(s); rm(ij)
  
  Delta_mat = Delta/no_pop
  Delta_mat_eigen = eigen(Delta_mat)
  
  # Two ways of selecting the number of components
  
  K_val_k_criterion = k_criterion(eigenvalue = Delta_mat_eigen$values, plot = FALSE)
  K_val_select_K = select_K(tau = 10^-3, eigenvalue = Delta_mat_eigen$values)
  
  # new way of selecting the number of components
  
  K_val = vector("numeric", sample_size)
  for(ik in 1:sample_size)
  {
    K_val[ik] = select_K_new(eigenvalue = Delta_mat_eigen$values, index = ik, sample_size = sample_size,
                             no_pop = no_pop)
    rm(ik)
  }
  
  q_val_est = min(max(c(which.min(K_val) - 1, K_val_k_criterion, K_val_select_K)),6)
  if(q_val_est == 0)
  {
    warning("The number of components is zero.")
  }
  Delta_eigen_vector = as.matrix(Delta_mat_eigen$vectors[,1:q_val_est]) * sqrt(sample_size)
  
  ###############################
  # estimate the eigen-dimension
  ###############################
  
  factor_loading = array(NA, dim = c(q_val_est, D_val, no_pop))
  for(ij in 1:no_pop)
  {
    factor_loading[,,ij] = (crossprod(Delta_eigen_vector, data[,ij,]))/sample_size
    rm(ij)
  }
  
  # variance
  
  C_X = array(NA, dim = c(no_pop, no_pop, D_val, D_val))
  for(ij in 1:no_pop)
  {
    for(iw in 1:no_pop)
    {
      C_X[ij,iw,,] = crossprod(factor_loading[,,ij], factor_loading[,,iw])
    }
  }
  rm(ij); rm(iw)
  
  error = array(NA, dim = c(sample_size, no_pop, D_val))
  for(ij in 1:no_pop)
  {
    error[,ij,] = data[,ij,] - (Delta_eigen_vector %*% factor_loading[,,ij])
    rm(ij)
  }
  return(list(raw_data = data, factor_loading = factor_loading,
              factors = Delta_eigen_vector, error = error))
}

