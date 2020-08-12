########################################################
# Demean data to eliminate fixed effects.
# Create Datamatrix for C-Lasso with time-fixed effects.
########################################################
library(pracma)

###########################################################################
# If only individual-fixed effects are present.
# Input:
#       X: regressor array of dimension N x T x p
#       y: dependent variable array of dimension N x p
# Output:
#       X_demean_vector: demeaned regressors in vector format (N*T x p)
#       Y_demean_vector: demeaned dependent variable in vector format (N*T x 1)
#       X_demean_original: demeaned regressors in original format (N x T x p)
#       Y_demean_original: demeaned dependent variable in original format (N x T x p)

indiv_demean = function(X, Y){
  N = dim(X)[1]
  T = dim(X)[2]
  p = dim(X)[3]

  # Eliminate individual fixed effects.
  X_bar = array(0, c(N, T, p))
  for (pp in 1:p) {
    X_bar[ , , pp] = repmat(apply(X, c(1,3), mean)[ , pp], 1, T)
  }
  X_demean = X - X_bar
  X_demean_vector = zeros(T*N, p)
  for (ppp in 1:p) {
    X_demean_vector[, ppp] = pracma::Reshape(t(X_demean[ , , ppp]), N*T, 1)
  }

  Y_demean = Y - repmat(array(rowMeans(Y), c(N,1)), 1, T)
  Y_demean_vector = pracma::Reshape(t(Y_demean), T*N, 1)

  return(list(X_demean_vector=X_demean_vector, 
              Y_demean_vector=Y_demean_vector,
              X_demean_original=X_demean,
              Y_demean_original=Y_demean))
}

##########################################################################
# If individual- and time-fixed effects are present.
# Input:
#       X: regressor array of dimension N x T x p
#       y: dependent variable array of dimension N x p
# Output:
#       X_demean_vector: demeaned regressors in vector format (N*T x p)
#       Y_demean_vector: demeaned dependent variable in vector format (N*T x 1)
#       X_demean_original: demeaned regressors in original format (N x T x p)
#       Y_demean_original: demeaned dependent variable in original format (N x T x p)

indiv_time_demean = function(X, Y){
  N = dim(X)[1]
  T = dim(X)[2]
  p = dim(X)[3]

  # Eliminate individual and time fixed effects.
  X_bar_time = array(0, c(N, T, p))
  X_bar_indiv = array(0, c(N, T, p))
  for (pp in 1:p) {
    X_bar_time[ , , pp] = repmat(apply(X, c(1,3), mean)[ , pp], 1, T)
    X_bar_indiv[ , , pp] = repmat(apply(X, c(2,3), mean)[ , pp], N, 1)
  }
  X_bar = array(0, c(N, T, p))
  for (pp in 1:p) {
    X_bar[, , pp] = repmat(apply(X_bar_indiv, c(1,3), mean)[ , pp], 1, T)
  }
  X_demean = X - X_bar_time - X_bar_indiv + X_bar
  X_demean_vector = zeros(T*N, p)     # NT x p
  for (ppp in 1:p) {
    X_demean_vector[, ppp] = pracma::Reshape(t(X_demean[ , , ppp]), N*T, 1)
  }

  Y_demean = Y - repmat(array(rowMeans(Y), c(N,1)), 1, T) - repmat(colMeans(Y), N, 1) + repmat(mean(rowMeans(Y)), N, T)
  Y_demean_vector = pracma::Reshape(t(Y_demean), T*N, 1)   # NT x 1

  return(list(X_demean_vector=X_demean_vector, 
              Y_demean_vector=Y_demean_vector,
              X_demean_original=X_demean,
              Y_demean_original=Y_demean))
}

############################################################################
# Stack regressors in matrix to write objective function in matrix format.
# Input:
#       X: regressor matrix of dimension N x T x p
# Output:
#       Datamatrix: regressor matrix used in C-Lasso with time-fixed effects (NT x Np)

datamatrix = function(X){
  N = dim(X)[1]
  T = dim(X)[2]
  p = dim(X)[3]

  X_bar = array(0, c(N, T, p))
  for (pp in 1:p) {
    X_bar[ , , pp] = repmat(apply(X, c(1,3), mean)[ , pp], 1, T)
  }
  X_demean = X - X_bar
  factors = kronecker(diag(N), ones(T,p)) - 1/N * ones(N*T, N*p)
  Xmatrix = array()
  for (i in 1:N) {
    Xi = pracma::Reshape(X_demean[i, , ], T, p)
    Xmatrix = cbind(Xmatrix, repmat(Xi, N, 1))
  }
  Xmatrix = Xmatrix[ , 2:(N*p+1)]
  Datamatrix = factors * Xmatrix         # NT x Np

  return(Datamatrix)
}
