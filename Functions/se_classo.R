#############################################################################
# Compute analytical standard errors of group-specific coefficient estimates.
# The formulas follow the asymptotic theory of Lu and Su (2017) and the 
# notation follows Su, Shi and Phillips (2013).
#############################################################################

library(pracma)
library(MASS)
library(classo)
library(matrixStats)

source("./Functions/calculate_residuals.R")
source("./Functions/demean_data.R")

##############################################################################
# Input:
#       X: regressor matrix of dimension N x T x p
#       Y: dependent variable of dimension N x T
#       K: number of groups
#       group: estimated groups for each individual (N x 1)
#       group_coef: estimated group coefficients (K x p)
#       time_effect: dummy indicating whether time-fixed effects included
# Output:
#       stde: standard errors of group-specific coefficient estimates (K x p)

std_error_est = function(X, Y, K, group, group_coef, time_effect) {
  N = dim(X)[1]
  T = dim(X)[2]
  p = dim(X)[3]
  
  stde = array(0, c(K, p))
  
  for(k in 1:K) {
    # Gather individuals within group k.
    group_k = (group == k)
    
    indiv = 1:N
    group_index = indiv[group_k]
    N_k = sum(group_k)
    
    # Demean data within group k.
    if(time_effect == TRUE) {
      if (N_k > 1) {
        X_k = X[group_index, , ]
        Y_k = Y[group_index, ]
        
        data_k_demean = indiv_time_demean(X=X_k, Y=Y_k)
        Y_k_demean_vector = data_k_demean$Y_demean_vector
        X_k_demean_vector = data_k_demean$X_demean_vector
        Datamatrix_k = datamatrix(X=X_k)
        
      } else {
        # Special case of one individual in group k.
        Y_k = array(Y[group_index, ], c(N_k, T))
        X_k = array(X[group_index, , ], c(N_k, T, p))
        
        data_k_demean = indiv_demean(X=X_k, Y=Y_k)
        Y_k_demean_vector = data_k_demean$Y_demean_vector
        X_k_demean_vector = data_k_demean$X_demean_vector 
        Datamatrix_k = X_k_demean_vector
      }
      
      # Calculate regression residuals within group k.
      uhat_k = calc_uhat_time(Datamatrix=Datamatrix_k, 
                              Y_demean_vector=Y_k_demean_vector, 
                              betai_hat=repmat(group_coef[k, ], N_k, 1), 
                              N=N_k, 
                              T=T, 
                              p=p)
    } else {
      # If no time-fixed effects are present.
      X_k = X[group_index, , ]
      Y_k = Y[group_index, ]
      
      data_k_demean = indiv_demean(X=X_k, Y=Y_k)
      X_k_demean_vector = data_k_demean$X_demean_vector
      Y_k_demean_vector = data_k_demean$Y_demean_vector 
      
      # Calculate regression residuals within group k.
      uhat_k = calc_uhat(X_demean_vector=X_k_demean_vector, 
                         Y_demean_vector=Y_k_demean_vector, 
                         betai_hat=repmat(group_coef[k, ], N_k, 1), 
                         N=N_k, 
                         T=T, 
                         p=p)
    }
    
    # Compute asymptotic variance and standard error of group estimate of group k.
    phi_k = 1/(N_k*T) * (t(X_k_demean_vector) %*% X_k_demean_vector)
    
    psi_k = zeros(p, p)
    for (i in 1:N_k) {
      for (t in 1:T) {
        for (s in 1:T) {
          psi_k = psi_k + (pracma::Reshape(X_k[i, t, ], p, 1) %*% t(pracma::Reshape(X_k[i, s, ], p, 1))) * uhat_k[i, t] * uhat_k[i, s]
        }
      }
    }
    psi_k = psi_k/(N_k*T)
    
    variance = ginv(phi_k) %*% psi_k %*% ginv(phi_k)
    stde[k, ] = t(sqrt(diag(variance)/(N_k*T)))
  }
  
  return(stde)
}