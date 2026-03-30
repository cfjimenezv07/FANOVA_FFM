###############################
# Numbers of factors selection
###############################

### The method of den Reijer et al. 2021 ###

# k_criterion <- function(eigenvalue, plot = TRUE)
# {
#     n = length(eigenvalue)
#     Hn = sum(1/(1:n))
# 
#     eigen_diff = eigenvalue[-n] - eigenvalue[-1]
#     lambda_bar = 1/((2:n)*Hn)
# 
#     n_cross = which(eigen_diff <= lambda_bar)[1]
# 
#     if(n_cross > 1)
#     {
#         n_return = n_cross - 1
#     } 
#     else
#     {
#         n_return = 1
#     }
# 
#     if(plot)
#     {
#         plot(eigen_diff[1:(n_cross+6)], type = "b", xlab = expression("Eigenvalue number" ~ italic(k)), ylab = expression("Eigenvalue" ~~~ lambda))
#         lines(lambda_bar[1:(n_cross+6)], type ="b", col = 2)
#         abline(v = n_return, lty = 2, col = 4)
#         legend("topright", lty = c(1,1), lwd = c(2,2), col = c(1,2), c(expression(lambda ~~ "(eigenvalue)  "), expression(bar(lambda) ~~ "(threshold)  ")))
#     }
#     return(n_return)
# }


k_criterion <-function(eigenvalue, plot = TRUE)
{
  n = length(eigenvalue)
  if (n < 2) return(1)
  
  Hn = sum(1/(1:n))
  eigen_diff = eigenvalue[-n] - eigenvalue[-1]
  lambda_bar = 1/((2:n)*Hn)
  
  # Find indices and handle the empty case to avoid NA
  temp_indices = which(eigen_diff <= lambda_bar)
  
  if(length(temp_indices) > 0) {
    n_cross = temp_indices[1]
  } else {
    n_cross = 1 # Fallback when no intersection is found
  }
  
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
    # Use min() to ensure we don't index beyond the available data length
    plot_limit = min(n - 1, n_cross + 6)
    
    plot(eigen_diff[1:plot_limit], type = "b", 
         xlab = expression("Eigenvalue number" ~ italic(k)), 
         ylab = expression("Eigenvalue" ~~~ lambda))
    lines(lambda_bar[1:plot_limit], type ="b", col = 2)
    abline(v = n_return, lty = 2, col = 4)
    legend("topright", lty = c(1,1), lwd = c(2,2), col = c(1,2), 
           c(expression(lambda ~~ "(eigenvalue)  "), expression(bar(lambda) ~~ "(threshold)  ")))
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
