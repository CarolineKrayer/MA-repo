############################################
# Bias correction with half-panel jackknife.
############################################
library(pracma)
library(MASS)
library(classo)
library(matrixStats)

source("./Functions/calculate_residuals.R")
source("./Functions/demean_data.R")

#####################################################
# Input:
#       X: regressor matrix (N x T x p)
#       Y: dependent variable (N x T)
#       N: cross-section dimension
#       T: time-series dimension
#       p: number of regressors
#       K: number of groups 
#       group: estimated group number for each individual (N x 1)
#       group_coef: estimated group coefficients (K x p)
#                   (can also be used for within estimator with K=1)
#       time_effect: dummy indicating whether time-fixed effects are included
# Output:
#       group_coef_corr: bias corrected estimates of group-specific parameters (K x p)

bias_correction = function(X, Y, N, T, p, K, group, group_coef, time_effect) {
  group_coef_corr = array(0, c(K, p))
  
  for (k in 1:K) {
    # Gather individuals within group k.
    group_k = (group == k)
    N_k = sum(group_k)
    
    if (N_k >= 2 * p/T) {
      indiv = 1:N
      group_index = indiv[group_k]
      
      # Demean data within group k.
      if(time_effect == TRUE) {
        if (N_k > 1) {
          Y_k = Y[group_index, ]
          X_k = X[group_index, , ]
          
          data_k_demean = indiv_time_demean(X=X_k, Y=Y_k)
          Y_k_demean_vector = data_k_demean$Y_demean_vector
          X_k_demean_vector = data_k_demean$X_demean_vector        
        } else {
          # Special case of one individual in group k.
          Y_k = array(Y[group_index, ], c(N_k, T))
          X_k = array(X[group_index, , ], c(N_k, T, p))
          
          data_k_demean = indiv_demean(X=X_k, Y=Y_k)
          Y_k_demean_vector = data_k_demean$Y_demean_vector
          X_k_demean_vector = data_k_demean$X_demean_vector  
        }
      } else {
        data_k_demean = indiv_demean(X=X_k, Y=Y_k)
        X_k_demean_vector = data_k_demean$X_demean_vector
        Y_k_demean_vector = data_k_demean$Y_demean_vector 
      }
      
      # Perform bias correction.
      bias_k = classo::SPJ_PLS(t=T, y=Y_k_demean_vector, x=X_k_demean_vector)
      a_k = lsfit(x=X_k_demean_vector, 
                  y=Y_k_demean_vector, 
                  intercept = FALSE)$coefficients
      group_coef_corr[k, ] = 2 * a_k - bias_k
      
    } else {
      group_coef_corr[k, ] = group_coef[k, ]
    }
  }
  
  return(group_coef_corr)
}

