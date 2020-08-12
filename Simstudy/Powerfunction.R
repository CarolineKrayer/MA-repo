##########################################################################
# Powerfunction of random coefficient model with growing disturbance term
##########################################################################
library(classo)
library(Rmosek)
library(pbapply)
library(parallel)
library(pracma)
library(MASS)
library(matrixStats)
library(Matrix)

rm(list=ls())
out_path = "./Simstudy/Results/"
source("./Functions/calculate_residuals.R")
source("./Functions/test_stat.R")
source("./Functions/DGP.R")
source("./Functions/demean_data.R")

# Count rejections of test.
# Input: 
#       data: generated regressors and dependent variable
#       alpha: significance level used for test 
#       time_effect: dummy indicating inclusion of time-fixed effects
#       bias_corr: dummy indicating whether bias correction should be performed
rejection_test = function(data, alpha, time_effect, bias_corr){
  X = data$X
  Y = data$Y
  N = dim(X)[1]
  T = dim(X)[2]
  p = dim(X)[3]
  
  # Demean and reshape data.
  if (time_effect == FALSE) {
    demean_data = indiv_demean(X=X, Y=Y)
    X_demean_vector = demean_data$X_demean_vector
    Y_demean_vector = demean_data$Y_demean_vector    
  } else {
    demean_data = indiv_time_demean(X=X, Y=Y)
    X_demean_vector = demean_data$X_demean_vector
    Y_demean_vector = demean_data$Y_demean_vector 
    Datamatrix = datamatrix(X=X)
  }

  # Compute common coefficient estimator.
  common_coef_est = lsfit(x=X_demean_vector, y=Y_demean_vector, intercept = FALSE)$coefficients
  common_coef_est = repmat(t(common_coef_est), N, 1)
  
  if (bias_corr == TRUE) {
    bias = SPJ_PLS(t=T, y=Y_demean_vector, x=X_demean_vector)
    common_coef_est = 2 * common_coef_est - repmat(bias, N, 1)
  } 
    
  # Compute regression residuals.
  if (time_effect == FALSE) {
    uhat_common = calc_uhat(X_demean_vector=X_demean_vector,
                            Y_demean_vector=Y_demean_vector,
                            betai_hat=common_coef_est,
                            N=N,
                            T=T,
                            p=p)
  } else {
    uhat_common = calc_uhat_time(Datamatrix=Datamatrix,
                                 Y_demean_vector=Y_demean_vector,
                                 betai_hat=common_coef_est,
                                 N=N,
                                 T=T,
                                 p=p)
  }
  
  # Compute test statistic for H0: K=1 and determine test decision.
  test_statistics = Jhat(uhat=uhat_common, X0=X)
  rej = as.numeric(test_statistics >= qnorm(1-alpha, 0, 1))
  
  return(rej)
}

###########################################################################
# Compute power of the test over Monte-Carlo repetitions.
# Input:
#       N: cross-sectional dimension
#       T: time series dimension
#       M: number of Monte Carlo repetitions
#       alpha: significance level used for test
#       time_effect: dummy indicating inclusion of time-fixed effects
#       bias_corr: dummy indicating whether bias correction should be performed
#       disturbance: standard deviation of disturbance term in coefficients

power_test = function(N, T, M, alpha, time_effect, bias_corr, disturbance) {
  # Generate data.
  set.seed(92642)
  data = lapply(rep(1, M), function(x) generate_data(DGP="power_DGP", 
                                                     time_effect=time_effect,
                                                     N=N, 
                                                     T=T, 
                                                     disturbance=disturbance))

  # Count rejections of the test.
  no_cores = detectCores() - 1 
  cl = makeCluster(no_cores)
  clusterEvalQ(cl, {
    library(classo)
    library(Rmosek)
    library(pracma)
    library(MASS)
    library(matrixStats)
    library(Matrix)
    library(expm)
  })
  clusterExport(cl=cl, 
                varlist=c("alpha", 
                          "bias_corr",
                          "indiv_demean", 
                          "indiv_time_demean",
                          "calc_uhat", 
                          "calc_uhat_time",
                          "Jhat",
                          "rejection_test",
                          "datamatrix"), 
                envir=environment())

  rejections = parSapply(cl, 
                         data, 
                         function(x) rejection_test(data=x, 
                                                    alpha=alpha, 
                                                    time_effect=time_effect,
                                                    bias_corr=bias_corr))
  stopCluster(cl)
  
  return(mean(rejections))
}

##################################################
# Simulation parameters.
N = 40
T = 10
M = 1000

time_effect = TRUE
bias_corr = TRUE
Kmax= 5
alpha = 0.05/Kmax
disturbances = seq(0, 0.4, by=0.005)

# Compute power of test.
power_results_dist = pbsapply(disturbances,
                              function(d) power_test(N=N,
                                                     T=T,
                                                     M=M,
                                                     alpha=alpha,
                                                     time_effect=time_effect,
                                                     bias_corr=bias_corr,
                                                     disturbance=d))

# Save results in a dataframe for later use. 
df_disturbances = data.frame(disturbances, power_results_dist)
names(df_disturbances) = c("disturbances", paste0("power_N", N, "_T", T))
saveRDS(df_disturbances, file=paste0(out_path, "Powerfct_homo_dist_N", N, "_T", T, "_Mrep", M, "_alpha", alpha, "_time_effect", time_effect, "_bias_corr", bias_corr, ".rds"))

