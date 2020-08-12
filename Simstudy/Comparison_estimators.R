#########################################################################
# Compute MSE, bias and empirical standard errors for various estimators.
#########################################################################
library(classo)
library(Rmosek)
library(pbapply)
library(parallel)
library(rlist)
library(pracma)
library(xtable)

rm(list=ls())
out_path = "./Simstudy/Results/"

#################################################################
# Compute MSE of estimator.
# Input:
#     beta_mat: matrix containing coefficient estimates (N x p)
#     beta_true: matrix with true coefficients (N x p)
# Output:
#     MSE: mean-squared error of coefficient estimates

calc_MSE = function(beta_mat, beta_true) {
  dev_beta_mat = mapply(function(x, y) x-y,
                        beta_mat,
                        beta_true,
                        SIMPLIFY=FALSE)
  sq_dev_beta_mat = lapply(dev_beta_mat, function(x) x**2)
  MSE = Reduce("+", sq_dev_beta_mat)/length(sq_dev_beta_mat)

  return(MSE)
}

#################################################################
# Compute bias of estimator.
# Input:
#     beta_mat: matrix containing coefficient estimates (N x p)
#     beta_true: matrix with true coefficients (N x p)
# Output:
#     bias: bias of coefficient estimates

calc_bias = function(beta_mat, beta_true) {
  dev_beta_mat = mapply(function(x, y) x-y,
                        x=beta_mat,
                        y=beta_true,
                        SIMPLIFY=FALSE)
  bias = Reduce("+", dev_beta_mat)/length(dev_beta_mat)
  
  return(bias)
}

#################################################################
# Compute empirical standard error of coefficient estimates.
# Input:
#     beta_mat: matrix containing estimated coefficients (N x p)
#     beta_true: matrix with true coefficients (N x p)
# Output:
#     emp_se: empirical standard error

calc_emp_se = function(beta_mat) {
  beta_mean = Reduce("+", beta_mat)/length(beta_mat)
  sq_dem_beta_mat = lapply(beta_mat, function(x) (x-beta_mean)**2)
  emp_var = Reduce("+", sq_dem_beta_mat)/(length(sq_dem_beta_mat)-1)
  emp_se = sqrt(emp_var)
  
  return(emp_se)
}


##########################################################################
# Gather results for comparison of MSE, bias and empirical standard error.
# Input:
#       N: cross-section dimension
#       T: time-series dimension
#       p: number of regressors
#       M: number of Monte-Carlo repetitions
#       Kmax: maximal number of groups considered
#       DGP: DGP considered
#       alpha: significance level used for test rejections
#       time_effect: dummy indicating whether time-fixed effects are included
#       ccons: constants for tuning parameter considered
#       bias_corr: dummy indicating whether bias correction should be performed
# Output:
#       MSE, bias and empirical standard error of estimators, averaged over
#       the coefficient estimates, respectively.

MSE_bias_results = function(N, T, p, M, Kmax, DGP, alpha, time_effect, ccons=c(0.25,0.5, 0.75), bias_corr=FALSE) {
  # Load results of simulation study.
  results = list.load(file=paste0(out_path, "Test_IC_results_DGP", DGP, "_time_effect", time_effect, "_N", N, "_T", T, "_Kmax", Kmax, "_Mrep", M, "_bias_corr", bias_corr,".rdata"))
  
  # Determine true beta.
  beta0 = results[7, ]
  
  # Common coefficient estimator.
  common_coef = results[3, ]
  MSE_common_coef = calc_MSE(beta_mat=common_coef, beta_true=beta0)
  bias_common_coef = calc_bias(beta_mat=common_coef, beta_true=beta0)
  emp_se_common_coef = calc_emp_se(beta_mat=common_coef)
  
  # Heterogeneous coefficient estimator.
  hetero_coef = results[4, ]
  MSE_hetero_coef = calc_MSE(beta_mat=hetero_coef, beta_true=beta0)
  bias_hetero_coef = calc_bias(beta_mat=hetero_coef, beta_true=beta0)
  emp_se_hetero_coef = calc_emp_se(beta_mat=hetero_coef)

  ## Post-Lasso estimator with number of groups determined by test.
  test_stats = results[1, ]
  info_criterions = results[2, ]
  post_lasso_coef = results[5, ]
  
  # Determine number of groups selected by test and IC, respectively.
  rej = lapply(test_stats, function(x) as.numeric(x >= qnorm(1-alpha, 0, 1)))
  rej_array = array(as.numeric(unlist(rej)), dim=c(Kmax, length(ccons), M))
  est_group_number = apply(rej_array, c(2, 3), function(x) min(which(x == 0)))  # length(ccons) x M
  info_criterion_array = array(as.numeric(unlist(info_criterions)), dim=c(length(ccons), M)) # length(ccons) x M
  
  # Initialise empty lists for results.
  beta_mat_test = vector(mode = "list", length = M)
  beta_mat_IC = vector(mode = "list", length = M)
  MSE_post_test_coef = vector(mode="list", length=length(ccons))
  MSE_post_IC_coef = vector(mode="list", length=length(ccons))
  bias_post_test_coef = vector(mode="list", length=length(ccons))
  bias_post_IC_coef = vector(mode="list", length=length(ccons))
  emp_se_post_test_coef = vector(mode="list", length=length(ccons))
  emp_se_post_IC_coef = vector(mode="list", length=length(ccons))
  for (j in 1:length(ccons)) {
    for(m in 1:M) {
      # Save estimates for number of groups selected by test and IC.
      if(is.infinite(est_group_number[j, m])) {
        beta_mat_test[[m]] = hetero_coef[[m]]
      } else if (est_group_number[j, m]==1) {
        beta_mat_test[[m]] = common_coef[[m]]
      } else {
        group_test = est_group_number[j, m]
        beta_mat_test[[m]] = post_lasso_coef[[m]][ , , group_test, j]
      }
      
      # Save results of information criterion.
      if(is.infinite(info_criterion_array[j, m])) {
        beta_mat_IC[[m]] = hetero_coef[[m]]
      } else if (info_criterion_array[j, m]==1) {
        beta_mat_IC[[m]] = common_coef[[m]]
      } else {
        group_IC = info_criterion_array[j, m]
        beta_mat_IC[[m]] = post_lasso_coef[[m]][ , , group_IC, j]
      }
    }
    
    # Results of testing procedure.
    MSE_post_test_coef[[j]] = calc_MSE(beta_mat=beta_mat_test, beta_true=beta0)
    bias_post_test_coef[[j]] = calc_bias(beta_mat=beta_mat_test, beta_true=beta0)
    emp_se_post_test_coef[[j]] = calc_emp_se(beta_mat=beta_mat_test)
    
    # Results for information criterion.
    MSE_post_IC_coef[[j]] = calc_MSE(beta_mat=beta_mat_IC, beta_true=beta0)
    bias_post_IC_coef[[j]] = calc_bias(beta_mat=beta_mat_IC, beta_true=beta0)
    emp_se_post_IC_coef[[j]] = calc_emp_se(beta_mat=beta_mat_IC)
  }

  # Infeasible oracle estimator.
  oracle_coef = results[6, ]
  MSE_oracle_coef = calc_MSE(beta_mat=oracle_coef, beta_true=beta0)
  bias_oracle_coef = calc_bias(beta_mat=oracle_coef, beta_true=beta0)
  emp_se_oracle_coef = calc_emp_se(beta_mat=oracle_coef)
  
  return(list(MSE_common_coef = MSE_common_coef,
              MSE_hetero_coef = MSE_hetero_coef,
              MSE_post_test_coef = MSE_post_test_coef,
              MSE_post_IC_coef = MSE_post_IC_coef,
              MSE_oracle_coef = MSE_oracle_coef,
              bias_common_coef = bias_common_coef,
              bias_hetero_coef = bias_hetero_coef,
              bias_post_test_coef = bias_post_test_coef,
              bias_post_IC_coef = bias_post_IC_coef,
              bias_oracle_coef = bias_oracle_coef,
              emp_se_common_coef = emp_se_common_coef,
              emp_se_hetero_coef = emp_se_hetero_coef,
              emp_se_post_test_coef = emp_se_post_test_coef,
              emp_se_post_IC_coef = emp_se_post_IC_coef,
              emp_se_oracle_coef = emp_se_oracle_coef)
         )
}

################################################################################
# Simulation parameters.
N_vec = c(40, 80, 40, 20)
T_vec = c(10, 10, 40, 40)
p = 2

ccons = c(0.25, 0.5, 0.75)
lambda_pos = 2

M = 1000
Kmax = 5
DGP = 6
time_effect = TRUE
bias_corr = TRUE
alpha = 0.05/Kmax

# Create table for comparison of estimators.
table = matrix(NA, nrow = 12, ncol = 5, dimnames = list(c("MSE N=40, T=10", "MSE N=80, T=10", "MSE N=40, T=40", "MSE N=20, T=40", "Bias N=40, T=10", "Bias N=80, T=10", "Bias N=40, T=40", "Bias N=20 , T=40", "Empir SE N=40, T=10", "Empir SE N=80, T=10", "Emp se N=40, T=40", "Empir SE N=20, T=40"), c("Within Estimator", "Heterogeneous Estimator", "C-Lasso with Test", "C-Lasso with IC", "Infeasible Estimator")))
  
for( i in 1:4){
  results = MSE_bias_results(N=N_vec[i], T=T_vec[i], p=p, M=M, Kmax=Kmax, DGP=DGP, alpha=alpha, ccons=ccons, bias_corr=bias_corr, time_effect = time_effect)
  
  table[i, 1] = mean(rowMeans(results$MSE_common_coef))
  table[i, 2] = mean(rowMeans(results$MSE_hetero_coef))
  table[i, 3] = mean(rowMeans(results$MSE_post_test_coef[[lambda_pos]]))
  table[i, 4] = mean(rowMeans(results$MSE_post_IC_coef[[lambda_pos]]))
  table[i, 5] = mean(rowMeans(results$MSE_oracle_coef))
  
  table[i+4, 1] = mean(rowMeans(results$bias_common_coef))
  table[i+4, 2] = mean(rowMeans(results$bias_hetero_coef))
  table[i+4, 3] = mean(rowMeans(results$bias_post_test_coef[[lambda_pos]]))
  table[i+4, 4] = mean(rowMeans(results$bias_post_IC_coef[[lambda_pos]]))
  table[i+4, 5] = mean(rowMeans(results$bias_oracle_coef))
  
  table[i+8, 1] = mean(rowMeans(results$emp_se_common_coef))
  table[i+8, 2] = mean(rowMeans(results$emp_se_hetero_coef))
  table[i+8, 3] = mean(rowMeans(results$emp_se_post_test_coef[[lambda_pos]]))
  table[i+8, 4] = mean(rowMeans(results$emp_se_post_IC_coef[[lambda_pos]]))
  table[i+8, 5] = mean(rowMeans(results$emp_se_oracle_coef))
}


table_tex = xtable(table, digits=3)
caption(table_tex) = paste0("Estimator Performance DGP", DGP, "time effect", time_effect, "biascorr", bias_corr, "\n")

print.xtable(x=table_tex,
             type="latex",
             file=paste0(out_path, "comparison_est_DGP", DGP, "_time_effect", time_effect, "_biascorr", bias_corr, "_Mrep", M, "_lambda_pos", lambda_pos,".tex"),
             caption.placement = "bottom",
             booktabs=TRUE
) 

