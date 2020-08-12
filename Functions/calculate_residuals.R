###################################
# Compute the regression residuals.
###################################
library(pracma)

# If only individual-fixed effects are present.
# Input:
#       X_demean_vector: demeaned regressor matrix in vector format (NT x p)
#       Y_demean_vector: demeaned dependent variable in vector format (NT x 1)
#       betai_hat: coefficient estimates (N x p)
#       N: cross section dimension
#       T: time dimension
#       p: number of regressors in model
# Output:
#       uhat: regression residuals (N x T)

calc_uhat = function(X_demean_vector, Y_demean_vector, betai_hat, N, T, p){

  indiv_dummy = kronecker(diag(N), ones(T, p))
  temp = repmat(X_demean_vector, 1, N) * indiv_dummy
  uhat_temp = Y_demean_vector - (temp %*% pracma::Reshape(t(betai_hat), N*p, 1))
  uhat = t(Reshape(uhat_temp, T, N))

  return(uhat)
}

# If time- and individual-fixed effects are present.
# Input:
#       Datamatrix: regressor matrix used in C-Lasso with time-fixed effects (NT x Np)
#       Y_demean_vector: demeaned dependent variable in vector format (NT x 1)
#       betai_hat: coefficient estimates (N x p)
#       N: cross section dimension
#       T: time dimension
#       p: number of regressors in model
# Output:
#       uhat: regression residuals (N x T)

calc_uhat_time = function(Datamatrix, Y_demean_vector, betai_hat, N, T, p){
  
  uhat_temp = Y_demean_vector - (Datamatrix %*% pracma::Reshape(t(betai_hat), N*p, 1))
  uhat = t(Reshape(uhat_temp, T, N))
  
  return(uhat)
}
