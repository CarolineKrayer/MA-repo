#############################
# Implement simulation study.
#############################
library(classo)
library(Rmosek)
library(pbapply)
library(parallel)
library(rlist)
library(pracma)
library(MASS)
library(matrixStats)
library(Matrix)
library(tidyverse)

rm(list=ls())
out_path = "./Simstudy/Results/"
source("./Functions/calculate_residuals.R")
source("./Functions/test_stat.R")
source("./Functions/DGP.R")
source("./Functions/demean_data.R")
source("./Functions/estimators.R")
source("./Functions/lasso_est_time.R")

# Input:
#       data: list containing regressors and dependent variable
#       DGP: DGP considered (integer between 1 and 8)
#       N: cross-section dimension
#       T: time-series dimension
#       Kmax: maximal number of groups considered
#       time_effect: dummy indicating inclusion of time-fixed effects
#       bias_corr: dummy indicating whether bias correction should be performed
#       ccons: constants for tuning parameter considered
# Output:
#       test_stat: test statistics computed (Kmax x length(ccons))
#       IC: number of groups determined by info criterion (length(ccons) x 1)
#       common_coef_est: common coefficient estimates (N x p)
#       heterogen_coef_est: heterogeneous coefficient estimates (N x p)
#       post_lasso_est: post-Lasso estimates (N x p x Kmax x length(ccons))
#       infeasible_est: infeasible oracle estimates (N x p)
#       beta0: matrix containing true coefficients (N x p)

simul_test = function(data, DGP, N, T, Kmax, time_effect, bias_corr=FALSE, ccons=c(0.25, 0.5, 0.75)){
  X = data$X
  Y = data$Y
  beta0 = data$beta
  p = dim(X)[3]
  
  # Initialise arrays for results.
  test_statistics = array(0, c(Kmax, length(ccons)))
  info_criterion = array(0, c(Kmax+1, length(ccons)))
  min_info_criterion = array(0, c(length(ccons), 1))
  common_coef_est = array(0, c(N, p))
  post_lasso_est = array(0, c(N, p, Kmax, length(ccons)))
  infeasible_est = array(0, c(N, p))

  # Demean and reshape data.
  if (time_effect == FALSE) {
    # Model without time fixed effects.
    demean_data = indiv_demean(X=X, Y=Y)
    X_demean_vector = demean_data$X_demean_vector
    Y_demean_vector = demean_data$Y_demean_vector
  } else {
    # Model with time fixed effects.
    demean_data = indiv_time_demean(X=X, Y=Y)
    X_demean_vector = demean_data$X_demean_vector
    Y_demean_vector = demean_data$Y_demean_vector
    Datamatrix = datamatrix(X=X)
  }

  ##################################################
  ################# Homogeneous model ##############
  ##################################################
  
  # Compute common coefficient estimator.
  common_coef_est = within_est(X_demean_vector=X_demean_vector,
                               Y_demean_vector=Y_demean_vector,
                               N=N,
                               T=T,
                               p=p,
                               bias_corr=bias_corr)$beta_hat
  
  # Compute regression residuals.
  if (time_effect == FALSE) {
    # Model without time fixed effects.
    uhat_common = calc_uhat(X_demean_vector=X_demean_vector,
                            Y_demean_vector=Y_demean_vector,
                            betai_hat=common_coef_est,
                            N=N,
                            T=T,
                            p=p)  
  } else {
    # Model with time fixed effects.
    uhat_common = calc_uhat_time(Datamatrix=Datamatrix,
                                 Y_demean_vector=Y_demean_vector,
                                 betai_hat=common_coef_est,
                                 N=N,
                                 T=T,
                                 p=p)
  }
  
  # Compute test statistic for H0: K=1.
  test_stat_common = Jhat(uhat=uhat_common, X0=X)
  test_statistics[1, ] = repmat(test_stat_common, 1, length(ccons))
  
  # Compute information criterion.
  sigma_hat_common = mean(rowMeans(uhat_common * uhat_common))
  rho = 2/3*(N*T)**(-1/2)
  info_criterion[1, ] = rep(log(sigma_hat_common) + rho*p, length(ccons))

  
  ##################################################
  ############# Group-heterogeneous model ##########
  ##################################################  
  
  lambda_vec = ccons*c(var(Y_demean_vector))*T**(-1/3)
  for (lambda_pos in 1:length(lambda_vec)) {
    if(lambda_pos==2) {
      K_stop = Kmax
    } else {
      K_stop = 3
    }
    for (k in 2:K_stop) {
      if (time_effect == FALSE) {
      # C-Lasso estimation for models without time fixed effects.
        beta_hat = PLS.mosek(N=N,
                             TT=T,
                             y=Y_demean_vector,
                             X=X_demean_vector,
                             K=k,
                             lambda=lambda_vec[lambda_pos],
                             bias_corr=bias_corr)$b.est
        post_lasso_est[ , , k, lambda_pos] = beta_hat
        uhat = calc_uhat(X_demean_vector=X_demean_vector,
                         Y_demean_vector=Y_demean_vector,
                         betai_hat=beta_hat,
                         N=N,
                         T=T,
                         p=p) 
      } else {
      # C-Lasso estimation for models with time fixed effects.
        beta_hat = lasso_est_time_mosek(X=X,
                                        Y=Y,
                                        X_demean_vector=X_demean_vector,
                                        Y_demean_vector=Y_demean_vector,
                                        Datamatrix=Datamatrix,
                                        K=k,
                                        lambda = lambda_vec[lambda_pos],
                                        bias_corr = bias_corr)$b.est
        post_lasso_est[ , , k, lambda_pos] = beta_hat
        uhat = calc_uhat_time(Datamatrix=Datamatrix,
                              Y_demean_vector=Y_demean_vector,
                              betai_hat=beta_hat,
                              N=N,
                              T=T,
                              p=p)
      }

      # Test H0: K = k.
      test_stat = Jhat(uhat=uhat, X0=X)
      test_statistics[k, lambda_pos] = test_stat

      # Compute information criterion.
      sigma_hat = mean(rowMeans(uhat * uhat))
      rho = 2/3*(N*T)**(-1/2)
      info_criterion[k, lambda_pos] = log(sigma_hat) + rho*p*k
    }
  }
  #####################################################
  ##### Compute other estimators for comparison. ######
  #####################################################
  # Heterogeneous coefficient estimator.
  heterogen_coef_est = hetero_coef_est(X=X, Y=Y)
  if (time_effect == FALSE) {
    uhat_hetero = calc_uhat(X_demean_vector=X_demean_vector,
                            Y_demean_vector=Y_demean_vector,
                            betai_hat=heterogen_coef_est,
                            N=N,
                            T=T,
                            p=p)
  } else {
    uhat_hetero = calc_uhat_time(Datamatrix=Datamatrix,
                                 Y_demean_vector=Y_demean_vector,
                                 betai_hat=heterogen_coef_est,
                                 N=N,
                                 T=T,
                                 p=p)
  }
  sigma_hat_hetero = mean(rowMeans(uhat_hetero * uhat_hetero))
  rho = 2/3*(N*(T-1))**(-1/2)
  info_criterion[Kmax+1, ] = rep(log(sigma_hat_hetero) + rho*p*N , length(ccons))
  
  # Infeasible oracle estimator.
  if (DGP == 7) {
    infeasible_est = heterogen_coef_est
  } else if (DGP == 8) {
    infeasible_est = swamy_est(X=X, Y=Y, N=N, T=T, p=p, time_effect=time_effect)
  } else {
    infeasible_est = oracle_est(X=X, Y=Y, DGP=DGP, time_effect=time_effect, bias_corr=bias_corr)
  }
  
  # Determine number of groups selected by information criterion.
  min_info_criterion = apply(info_criterion, 2, function(x) which.min(x))

  return(list(test_stat=test_statistics,
              IC=min_info_criterion,
              common_coef_est=common_coef_est,
              heterogen_coef_est=heterogen_coef_est,
              post_lasso_est=post_lasso_est,
              infeasible_est=infeasible_est,
              beta0=beta0))
}

################################################################################
# Run the simulation study.
N = 40
T = 40
M = 1000
Kmax = 5
DGP = 6
ccons=c(0.25, 0.5, 0.75)
bias_corr = TRUE
time_effect = TRUE

start = Sys.time()
# Generate data.
set.seed(92642)
data = lapply(rep(1, M), 
              function(x) generate_data(DGP=DGP, 
                                        time_effect=time_effect, 
                                        N=N, 
                                        T=T))

# Parallelise simulation study.
no_cores = detectCores() - 1 
cl = makeCluster(no_cores)
clusterEvalQ(cl, {
  library(classo)
  library(Rmosek)
  library(rlist)
  library(pracma)
  library(MASS)
  library(matrixStats)
  library(Matrix)
})
clusterExport(cl=cl, 
              varlist=c("simul_test", 
                        "Kmax", 
                        "N", 
                        "T", 
                        "M", 
                        "DGP", 
                        "bias_corr",
                        "time_effect",
                        "ccons",
                        "generate_data", 
                        "indiv_demean", 
                        "indiv_time_demean",
                        "datamatrix",
                        "within_est",
                        "calc_uhat", 
                        "calc_uhat_time",
                        "Jhat", 
                        "hetero_coef_est", 
                        "oracle_est",
                        "lasso_est_time_mosek",
                        "pen_generate",
                        "opt_mosek_time",
                        "criterion_time",
                        "post_corr_time"), 
              envir=environment())

output = parSapply(cl, 
                   data, 
                   function(x) simul_test(data = x,
                                          DGP = DGP,
                                          N = N,
                                          T = T,
                                          Kmax = Kmax,
                                          time_effect = time_effect,
                                          bias_corr = bias_corr,
                                          ccons = ccons))

stopCluster(cl)
end = Sys.time()
runtime = end - start

# Save results in a list.
list.save(output, file=paste0(out_path, "Test_IC_results_DGP", DGP, "_time_effect", time_effect, "_N", N, "_T", T, "_Kmax", Kmax, "_Mrep", M, "_bias_corr", bias_corr,".rdata"))
