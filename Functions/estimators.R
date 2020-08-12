#########################################################
# Implementation of estimators apart from C-Lasso
#########################################################
library(pracma)
library(plm)
library(MASS)
library(classo)
library(stats)

source("./Functions/demean_data.R")
source("./Functions/calculate_residuals.R")

###############################################################################
# Fixed-effects estimator.
# Input:
#       X_demean_vector: demeaned regressor matrix in vector format (NT x p)
#       Y_demean_vector: demeaned dependent variable in vector format (NT x 1)
#       N: cross section dimension
#       T: time dimension
#       p: number of regressors in model
#       bias_corr: dummy whether bias correction should be performed
# Output:
#       beta_hat: matrix with (bias corrected) coefficient estimates (N x p)
#       se: standard error of coefficient estimates in matrix format (N x p)

within_est = function(X_demean_vector, Y_demean_vector, N, T, p, bias_corr=FALSE){
  reg = lsfit(x=X_demean_vector, y=Y_demean_vector, intercept = FALSE)
  beta_hat = repmat(t(reg$coefficients), N, 1)
  se = ls.diag(reg)$std.err
  
  if (bias_corr == TRUE && N >= 2*p/T) {
    # Bias correction with half-panel jackknife.
    bias = SPJ_PLS(t=T, y=Y_demean_vector, x=X_demean_vector)
    beta_hat_corr = 2 * beta_hat - repmat(bias, N, 1)
    
    return(list(beta_hat = beta_hat_corr,
                se = se))
  } else {
    return(list(beta_hat = beta_hat,
                se = se))
  }
}

###########################################################################
# Completely heterogeneous coefficient estimator.
# Input:
#       X: regressor array of dimension N x T x p
#       Y: dependent variable array of dimension N x T
# Output:
#       beta_hetero_coef: estimated coefficients (N x p)


hetero_coef_est = function(X, Y){
  N = dim(X)[1]
  T = dim(X)[2]
  p = dim(X)[3]
  
  beta_hetero_coef = zeros(N, p)

  for (i in 1:N) {
    regressors = cbind(pracma::Reshape(X[i, , ], T, p), ones(T, 1))
    dep_var = Y[i, ]
    reg = lsfit(x=regressors, y=dep_var, intercept=FALSE)
    beta_est = reg$coefficients

    beta_hetero_coef[i, ] = t(beta_est[1:p])      # N x p
  }
  
  return(beta_hetero_coef)
}

###########################################################################
# Infeasible oracle estimator.
# Input:
#       X: regressor array of dimension N x T x p
#       Y: dependent variable array of dimension N x T
#       DGP: DGP used to generate data
#       bias_corr: dummy indicating whether bias correction should be performed
# Output:
#       infeasible beta/ infeasible_beta_corr: (bias corrected) coefficient estimates (N x p)

oracle_est = function(X, Y, DGP, time_effect, bias_corr=FALSE){
  N = dim(X)[1]
  T = dim(X)[2]
  p = dim(X)[3]
  
  infeasible_beta = array()
  infeasible_beta_corr = array()
  
  if (DGP != 7 && DGP != 8) {
    # Determine number of individuals in groups.
    if (DGP == 5){
      N1 = floor(N*0.7)
      N2 = floor(N*0.05)
    } else{
      N1 = floor(N*0.3)
      N2 = floor(N*0.3)
    }

    # Minimal number of groups needed for bias correction.
    group_size_criterion = 0
    
    # Estimate coefficients within the true groups separately. 
    for (group in 1:3) {
      # Determine group membership.
      if (group == 1) {
        group_X_true = array(X[1:N1, , ], c(N1, T, p))
        group_Y_true = array(Y[1:N1, ], c(N1, T))
      }
      if (group == 2) {
        group_X_true = array(X[(N1+1):(N1+N2), , ], c(N2, T, p))
        group_Y_true = array(Y[(N1+1):(N1+N2), ], c(N2, T))
      }
      if (group == 3) {
        group_X_true = array(X[(N1+N2+1):N, , ], c(N-N1-N2, T, p))
        group_Y_true = array(Y[(N1+N2+1):N, ], c(N-N1-N2, T))
      }
      
      # Demean data within group.
      if (time_effect == FALSE) {
        demean_group_data = indiv_demean(X=group_X_true, Y=group_Y_true)
        group_X_true_demean_vector = demean_group_data$X_demean_vector
        group_Y_true_demean_vector = demean_group_data$Y_demean_vector
      } else {
        demean_group_data = indiv_time_demean(X=group_X_true, Y=group_Y_true)
        group_X_true_demean_vector = demean_group_data$X_demean_vector
        group_Y_true_demean_vector = demean_group_data$Y_demean_vector
      }
      
      # Estimate coefficients within group.
      N_group = dim(group_X_true)[1]
      beta_hat = within_est(X_demean_vector=group_X_true_demean_vector,
                            Y_demean_vector=group_Y_true_demean_vector,
                            N=N_group,
                            T=T,
                            p=p,
                            bias_corr=FALSE)$beta_hat
      infeasible_beta = rbind(infeasible_beta, beta_hat)
      
      if (bias_corr==TRUE && N_group >= 2*p/T) {
        # Bias correction with half-panel jackknife.
        group_size_criterion = group_size_criterion + 1
        bias_group = SPJ_PLS(t=T,
                             y=group_Y_true_demean_vector,
                             x=group_X_true_demean_vector)
        beta_hat_corr = 2 * beta_hat - repmat(bias_group, N_group, 1)
        
        infeasible_beta_corr = rbind(infeasible_beta_corr, beta_hat_corr)
      }
    }
  } else {
    print("No group structure in coefficients - reconsider implementation of DGP.")
  }
  
  
  if (bias_corr==TRUE && group_size_criterion==3) {
    infeasible_beta_corr = infeasible_beta_corr[2:(N+1), ]
    return(infeasible_beta_corr)
  } else {
    infeasible_beta = infeasible_beta[2:(N+1), ]
    return(infeasible_beta)
  }
}

#################################################################################
# Swamy's random coefficient estimator -formulas follow Hsiao (2014, pp. 172-175)
# Input:
#       X: regressor matrix of dimension N x T x p
#       Y: dependent variable of dimension N x T
#       N: cross-sectional dimension
#       T: time series dimension
#       p: number of regressors
#       time_effect: dummy indicating whether model includes time effects
# Output:
#       swamy_coefs: Swamy's random coefficient estimator (p x 1)
#       indiv_swamy_coefs: individual estimates computed from Swamy's estimator (N x p)

swamy_est = function(X, Y, N, T, p, time_effect) {
  # Demean and reshape X and Y.
  if (time_effect == TRUE) {
    demean_data = indiv_time_demean(X=X, Y=Y)
    X_demean_vector = demean_data$X_demean_vector
    Y_demean_vector = demean_data$Y_demean_vector
    X_demean_original = demean_data$X_demean_original
    Y_demean_original = demean_data$Y_demean_original
  } else {
    demean_data = indiv_demean(X=X, Y=Y)
    X_demean_vector = demean_data$X_demean_vector
    Y_demean_vector = demean_data$Y_demean_vector 
    X_demean_original = demean_data$X_demean_original
    Y_demean_original = demean_data$Y_demean_original
  }
  
  # Compute completely heterogeneous coefficient estimates. 
  beta_i = hetero_coef_est(X=X, Y=Y)
  beta_i_res = calc_uhat(X_demean_vector=X_demean_vector, 
                         Y_demean_vector=Y_demean_vector, 
                         betai_hat=beta_i, 
                         N=N, 
                         T=T, 
                         p=p)
  
  # Estimate variance of error terms.
  sigma_i_hat = diag(beta_i_res %*% t(beta_i_res))/(T-p)
  
  # Estimate variance of disturbance term in coefficients.
  first_term = array(0, c(p, p))
  sec_term = array(0, c(p, p))
  for (i in 1:N) {
    X_i = pracma::Reshape(X_demean_original[i, , ], T, p)
    first_term = first_term + ((beta_i[i, ] - 1/N * colSums(beta_i)) %*% t(beta_i[i, ] - 1/N * colSums(beta_i)))
    sec_term = sec_term + (sigma_i_hat[i] * ginv(t(X_i) %*% X_i))
  }
  
  first_term = first_term/(N - 1)
  sec_term = sec_term/N
  delta_hat = first_term - sec_term
  
  if(min(eigen(delta_hat)$values)<0) {
    delta_hat = first_term
  }
  
  # Compute Swamy's random coefficient estimate.
  temp1 = array(0, c(p, p))
  for (i in 1:N) {
    X_i = pracma::Reshape(X_demean_original[i, , ], T, p)
    temp1 = temp1 + ginv(delta_hat + sigma_i_hat[i] * ginv(t(X_i) %*% X_i))
  }
  
  swamy_coefs = array(0, c(p, 1))
  for (j in 1:N) {
    X_j = pracma::Reshape(X[j, , ], T, p)
    swamy_coefs = swamy_coefs + ginv(temp1) %*% ginv(delta_hat + sigma_i_hat[j] * ginv(t(X_j) %*% X_j)) %*% beta_i[j, ]
  }
  
  # Compute individual specific coefficient estimates.
  indiv_swamy_coefs = array(0, c(N, p))
  for(i in 1:N) {
    X_i = pracma::Reshape(X_demean_original[i, , ], T, p)
    y_i = pracma::Reshape(Y_demean_original[i, ], T, 1)
    
    indiv_swamy_coefs[i, ] = swamy_coefs + delta_hat %*% t(X_i) %*% ginv(X_i %*% delta_hat %*% t(X_i) + sigma_i_hat[i] * eye(T)) %*% (y_i - X_i %*% swamy_coefs)
  }
  
  return(list(swamy_coefs=swamy_coefs,
              indiv_swamy_coefs=indiv_swamy_coefs))
}