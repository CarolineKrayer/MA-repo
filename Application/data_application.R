#########################################################
# Function that computes estimates for data application.
#########################################################
library(pracma)
library(dplyr)
library(devtools)
library(plm)
library(classo)
library(Rmosek)
library(MASS)
library(matrixStats)
library(Matrix)
library(tidyverse)
library(rlist)
library(tis)

source("./Functions/calculate_uhat.R")
source("./Functions/test_stat.R")
source("./Functions/demean_data.R")
source("./Functions/estimators.R")
source("./Functions/lasso_est_time.R")
source("./Functions/bias_variance.R")

################################################################################
# Input:
#       X: regressor matrix of dimension N x T x p
#       Y: dependent variable of dimension N x T
#       N: cross-sectional dimension
#       T: time series dimension
#       p: number of regressors
#       Kmax: maximal number of groups considered
#       time_effect: dummy indicating inclusion of time-fixed effects 
#                    in model
#       bias_corr: dummy indicating whether bias correction is performed
#       ccons: constants for tuning parameter lambda considered
#
# Output:
#       test_stat: test statistics computed (Kmax x length(ccons))
#       p_test_stat: p-value of test statistic (Kmax x length(ccons))
#       min_IC: number of groups determined by info criterion 
#               (length(ccons) x 1)
#       IC: value of information criterion (Kmax+1 x length(ccons))
#       common_coef_est: common coefficient estimator (N x p)
#       common_coef_est_corr: bias corrected common coefficient estimator (N x p)
#       post_lasso_est: post-Lasso estimators (Kmax x p x Kmax x length(ccons))
#       post_lasso_est_corr: bias corrected post-Lasoo estimators 
#                            (Kmax x p x Kmax x length(ccons))
#       group: group membership of each individual estimated with C-Lasso 
#              (N x Kmax x length(ccons))
#       se: standard error of coefficient estimates (Kmax x p x Kmax x length(ccons))
#       t_test: test for significance of coefficient estimates (Kmax x p x Kmax x length(ccons))
#       p_estim: p-value of t-test (Kmax x p x Kmax x length(ccons))
#       hetero_coef_est: heterogeneous coefficient estimates (N x p)
            
data_application = function(X, Y, N, T, p, Kmax, time_effect, bias_corr, ccons=c(0.25,0.5,0.75)) {
  # Demean data.
  if(time_effect == TRUE) {
    data_demean = indiv_time_demean(X=X, Y=Y)
    X_demean_vector = data_demean$X_demean_vector
    Y_demean_vector = data_demean$Y_demean_vector
    Datamatrix = datamatrix(X=X)    
  } else {
    data_demean = indiv_demean(X=X, Y=Y)
    X_demean_vector = data_demean$X_demean_vector
    Y_demean_vector = data_demean$Y_demean_vector   
  }
  
  # Initialise arrays for results.
  test_statistics = zeros(Kmax, length(ccons))
  p_test_stat = zeros(Kmax, length(ccons))
  info_criterion = array(0, c(Kmax+1, length(ccons)))
  min_info_criterion = array(0, c(length(ccons), 1))
  post_lasso_est = array(0, c(Kmax, p, Kmax, length(ccons)))
  post_lasso_est_corr = array(0, c(Kmax, p, Kmax, length(ccons)))
  group = array(0, c(N, Kmax, length(ccons)))
  se = array(0, c(Kmax, p, Kmax, length(ccons)))
  t_test = array(0, c(Kmax, p, Kmax, length(ccons)))
  p_estim = array(0, c(Kmax, p, Kmax, length(ccons)))
  
  # Fix parameters needed for C-Lasso.
  rho = 2/3*(N*T)**(-1/2)
  beta0 = matrix(rnorm(N*p), N, p)
  
  ####################
  # Homogeneous model. 
  ####################
  # Within estimation.
  within_reg = within_est(X_demean_vector=X_demean_vector, 
                          Y_demean_vector=Y_demean_vector, 
                          N=N, 
                          T=T, 
                          p=p)
  
  common_coef_est = within_reg$beta_hat
  se[1, , 1, ] = repmat(within_reg$se, 1, length(ccons))
  
  if (time_effect == TRUE) {
    # Compute residuals in model with time-fixed effects.
    uhat_common_wo_bias_corr = calc_uhat_time(Datamatrix=Datamatrix,
                                              Y_demean_vector=Y_demean_vector,
                                              betai_hat=common_coef_est,
                                              N=N,
                                              T=T,
                                              p=p)
    
    if (bias_corr == TRUE) {
      # Bias correction with half-panel jackknife.
      common_coef_est_corr = bias_correction(X=X, 
                                             Y=Y, 
                                             N=N, 
                                             T=T, 
                                             p=p, 
                                             K=1, 
                                             group=c(rep(1, N)), 
                                             group_coef=common_coef_est,
                                             time_effect=time_effect)
      
      uhat_common = calc_uhat_time(Datamatrix=Datamatrix,
                                   Y_demean_vector=Y_demean_vector,
                                   betai_hat=repmat(common_coef_est_corr, N, 1),
                                   N=N,
                                   T=T,
                                   p=p)
      
      # Test significance of coefficient estimates.
      t_test[1, , 1, ] = repmat(t(common_coef_est_corr), 1, length(ccons))/se[1, , 1, ]
      p_estim[1, , 1, ] = 1-pnorm(t_test[1, , 1, ])
      
    } else {
      uhat_common = uhat_common_wo_bias_corr
      
      # Test significance of coefficient estimates.
      t_test[1, , 1, ] = t(repmat(common_coef_est[1, ], length(ccons), 1))/se[1, , 1, ]
      p_estim[1, , 1, ] = 1-pnorm(t_test[1, , 1, ])
    }
    
  } else {
    # Compute residuals in model without time-fixed effects.
    uhat_common_wo_bias_corr = calc_uhat(X_demean_vector=X_demean_vector, 
                                         Y_demean_vector=Y_demean_vector, 
                                         betai_hat=common_coef_est, 
                                         N=N, 
                                         T=T, 
                                         p=p)
    
    if (bias_corr == TRUE) {
      # Bias correction with half-panel jackknife.
      common_coef_est_corr = bias_correction(X=X, 
                                             Y=Y, 
                                             N=N, 
                                             T=T, 
                                             p=p, 
                                             K=1, 
                                             group=c(rep(1, N)), 
                                             group_coef=common_coef_est,
                                             time_effect=time_effect)
      
      uhat_common = calc_uhat(X_demean_vector=X_demean_vector, 
                              Y_demean_vector=Y_demean_vector, 
                              betai_hat=repmat(common_coef_est_corr, N, 1), 
                              N=N, 
                              T=T, 
                              p=p)
      
      # Test significance of coefficient estimates.
      t_test[1, , 1, ] = repmat(t(common_coef_est_corr), 1, length(ccons))/se[1, , 1, ]
      p_estim[1, , 1, ] = 1-pnorm(t_test[1, , 1, ])
      
    } else {
      uhat_common = uhat_common_wo_bias_corr
      
      # Test significance of coefficient estimates.
      t_test[1, , 1, ] = t(repmat(common_coef_est[1, ], length(ccons), 1))/se[1, , 1, ]
      p_estim[1, , 1, ] = 1-pnorm(t_test[1, , 1, ])
    }
  }
  
  # Test H0: K0=1.
  test_stat_common = Jhat(uhat=uhat_common, X0=X)
  test_statistics[1, ] = repmat(test_stat_common, 1, length(ccons))
  p_test_stat [1, ] = rep(1-pnorm(test_stat_common), length(ccons))
  
  # Compute information criterion.
  sigma_hat_common = mean(rowMeans(uhat_common * uhat_common))
  info_criterion[1, ] = rep(log(sigma_hat_common) + rho*p, length(ccons))
  
  #################################################
  # Group heterogeneous model - C-Lasso estimation.
  #################################################
  lambda_vec = ccons*c(var(Y_demean_vector))*T**(-1/3)
  
  for (lambda_pos in 1:length(lambda_vec)) {
    for (k in 2:Kmax) {
      if (time_effect == TRUE) {
        # C-Lasso estimation with time-fixed effects.
        lasso_est = lasso_est_time_mosek(X=X,
                                         Y=Y,
                                         X_demean_vector=X_demean_vector,
                                         Y_demean_vector=Y_demean_vector,
                                         Datamatrix=Datamatrix,
                                         K=k,
                                         lambda=lambda_vec[lambda_pos],
                                         beta0=beta0)
        
        beta_hat = lasso_est$b.est
        group_coef = lasso_est$a.out
        group[ , k, lambda_pos] = lasso_est$group.est
        post_lasso_est[1:k, , k, lambda_pos] = group_coef
        
        uhat_wo_bias_corr = calc_uhat_time(Datamatrix=Datamatrix,
                                           Y_demean_vector=Y_demean_vector,
                                           betai_hat=beta_hat,
                                           N=N,
                                           T=T,
                                           p=p)
        
        se[1:k, , k, lambda_pos] = std_error_est(X=X, 
                                                 Y=Y, 
                                                 K=k, 
                                                 group=group[ , k, lambda_pos], 
                                                 group_coef=group_coef, 
                                                 time_effect=time_effect) 
        
        if (bias_corr == TRUE) {
          # Bias correction with half-panel jackknife.
          group_coef_corr = bias_correction(X=X, 
                                            Y=Y, 
                                            N=N, 
                                            T=T, 
                                            p=p, 
                                            K=k, 
                                            group=group[ , k, lambda_pos], 
                                            group_coef=group_coef,
                                            time_effect=time_effect)
          post_lasso_est_corr[1:k, , k, lambda_pos] = group_coef_corr
          
          # Determine individual coefficients.
          beta_hat_corr = zeros(N, p)
          for (i in 1:N) {
            group_i = group[i, k, lambda_pos]
            beta_hat_corr[i, ] = post_lasso_est_corr[group_i, , k, lambda_pos]
          }
          
          uhat = calc_uhat_time(Datamatrix=Datamatrix,
                                Y_demean_vector=Y_demean_vector,
                                betai_hat=beta_hat_corr,
                                N=N,
                                T=T,
                                p=p)
          
          # Test significance of coefficient estimates.
          t_test[1:k, , k, lambda_pos] = post_lasso_est_corr[1:k, , k, lambda_pos]/se[1:k, , k, lambda_pos]
          p_estim[1:k, , k, lambda_pos] = 1-pnorm(t_test[1:k, , k, lambda_pos])
          
        } else {
          uhat = uhat_wo_bias_corr
          
          # Test significance of coefficient estimates.
          t_test[1:k, , k, lambda_pos] = post_lasso_est[1:k, , k, lambda_pos]/se[1:k, , k, lambda_pos]
          p_estim[1:k, , k, lambda_pos] = 1-pnorm(t_test[1:k, , k, lambda_pos])
        }
        
      } else {
        # C-Lasso estimation without time-fixed effects. 
        lasso_est = PLS.mosek(N=N,
                              TT=T,
                              y=Y_demean_vector,
                              X=X_demean_vector,
                              K=k,
                              lambda=lambda_vec[lambda_pos],
                              beta0=beta0)
        
        beta_hat = lasso_est$b.est
        group_coef = lasso_est$a.out
        group[ , k, lambda_pos] = lasso_est$group.est
        post_lasso_est[1:k, , k, lambda_pos] = group_coef
        
        uhat_wo_bias_corr = calc_uhat(X_demean_vector=X_demean_vector,
                                      Y_demean_vector=Y_demean_vector,
                                      betai_hat=beta_hat,
                                      N=N,
                                      T=T,
                                      p=p)
        
        se[1:k, , k, lambda_pos] = std_error_est(X=X, 
                                                 Y=Y, 
                                                 K=k, 
                                                 group=group[ , k, lambda_pos], 
                                                 group_coef=group_coef, 
                                                 time_effect=time_effect) 
        
        if (bias_corr == TRUE) {
          # Bias correction with half-panel jackknife.
          group_coef_corr = bias_correction(X=X, 
                                            Y=Y, 
                                            N=N, 
                                            T=T, 
                                            p=p, 
                                            K=k, 
                                            group=group[ , k, lambda_pos], 
                                            group_coef=group_coef,
                                            time_effect=time_effect)
          post_lasso_est_corr[1:k, , k, lambda_pos] = group_coef_corr
          
          # Determine individual coefficients.
          beta_hat_corr = zeros(N, p)
          for (i in 1:N) {
            group_i = group[i, k, lambda_pos]
            beta_hat_corr[i, ] = post_lasso_est_corr[group_i, , k, lambda_pos]
          }
          
          uhat = calc_uhat(X_demean_vector=X_demean_vector,
                           Y_demean_vector=Y_demean_vector,
                           betai_hat=beta_hat_corr,
                           N=N,
                           T=T,
                           p=p)
          
          # Test significance of coefficient estimates.
          t_test[1:k, , k, lambda_pos] = post_lasso_est_corr[1:k, , k, lambda_pos]/se[1:k, , k, lambda_pos]
          p_estim[1:k, , k, lambda_pos] = 1-pnorm(t_test[1:k, , k, lambda_pos])
          
        } else {
          uhat = uhat_wo_bias_corr
          
          # Test significance of coefficient estimates.
          t_test[1:k, , k, lambda_pos] = post_lasso_est[1:k, , k, lambda_pos]/se[1:k, , k, lambda_pos]
          p_estim[1:k, , k, lambda_pos] = 1-pnorm(t_test[1:k, , k, lambda_pos])
        }
      }
      
      # Test H0: K0 = k.
      test_stat = Jhat(uhat=uhat, X0=X)
      test_statistics[k, lambda_pos] = test_stat
      p_test_stat[k, lambda_pos] = 1-pnorm(test_stat)
      
      # Compute information criterion.
      sigma_hat = mean(rowMeans(uhat * uhat))
      info_criterion[k, lambda_pos] = log(sigma_hat) + rho*p*k
    }
  }

  #####################################
  # Completely heterogeneous estimator.
  #####################################
  heterogen_coef_est = hetero_coef_est(X=X, Y=Y)
  
  if (time_effect==TRUE) {
    uhat_hetero = calc_uhat_time(Datamatrix=Datamatrix,
                                 Y_demean_vector=Y_demean_vector,
                                 betai_hat=heterogen_coef_est,
                                 N=N,
                                 T=T,
                                 p=p)
  } else {
    uhat_hetero = calc_uhat(X_demean_vector=X_demean_vector,
                            Y_demean_vector=Y_demean_vector,
                            betai_hat=heterogen_coef_est,
                            N=N,
                            T=T,
                            p=p)
  }
  sigma_hat_hetero = mean(rowMeans(uhat_hetero * uhat_hetero))
  info_criterion[Kmax+1, ] = rep(log(sigma_hat_hetero) + rho*p*N, length(ccons))
  
  
  # Find K with smallest information criterion.
  for (lambda_pos in 1:length(lambda_vec)) {
    min_info_criterion[lambda_pos, 1] = which.min(info_criterion[ , lambda_pos])
  }
  
  # Return results.
  if(bias_corr == TRUE) {
    return(list(test_stat = test_statistics,
                p_test_stat = p_test_stat,
                min_IC = min_info_criterion,
                IC = info_criterion,
                common_coef_est = common_coef_est,
                common_coef_est_corr = common_coef_est_corr,
                post_lasso_est = post_lasso_est,
                post_lasso_est_corr = post_lasso_est_corr,
                group = group,
                se = se,
                t_test = t_test,
                p_estim = p_estim,
                hetero_coef_est = heterogen_coef_est))    
  } else {
    return(list(test_stat = test_statistics,
                p_test_stat = p_test_stat,
                min_IC = min_info_criterion,
                IC = info_criterion,
                common_coef_est = common_coef_est,
                post_lasso_est = post_lasso_est,
                group = group,
                se = se,
                t_test = t_test,
                p_estim = p_estim, 
                hetero_coef_est=heterogen_coef_est)) 
  }
}

