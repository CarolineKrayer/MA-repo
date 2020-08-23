################################################
# Implement hypothesis test of Lu and Su (2017).
################################################
library(MASS)
library(matrixStats)
library(pracma)
library(expm)

# Input:
#       uhat: matrix containing regression residuals (N x T)
#       X0: regressor array of dimension N x T x p
# Output:
#       Jhat: standardised LM-type test statistic (scalar)

Jhat = function(uhat, X0){
  N = dim(X0)[1]
  T = dim(X0)[2]
  p = dim(X0)[3]

  X = array(0, c(p,T,N))      # p x T x N
  for (ii in 1:p) {
    X[ii, , ] = t(X0[, , ii])
  }

  M0 = diag(T) - 1/T * array(1, c(T,T))    # T x T

  # Check whether any Omegaj_hat is singular and (if yes) delete individual j.
  checkproblems = array(0, c(N,1))
  for (j in 1:N) {
    Xj = t(X[ , ,j])
    Omegaj_hat = t(Xj) %*% M0 %*% Xj/T
    checkproblems[j,] = as.numeric(min(eigen(Omegaj_hat)$values) < 10**(-4))
  }
  if (length(which(checkproblems==1))!=0) {
    X = X[, , -which(checkproblems==1)]
    uhat = uhat[-which(checkproblems==1), ]
  }
  N = N - sum(checkproblems)

  # Calculate test statistic
  LM0 = array(0, c(N,1))
  B0 = array(0, c(N,1))
  V0 = array(0, c(N,1))

  for (i in 1:N) {
    # Raw statistic.
    Xi = t(X[, ,i])     # T x p
    if(p==1){
      Xi = t(Xi)
    }
    Omegai_hat = t(Xi) %*% M0 %*% Xi/T      # p x p
    PXi = M0 %*% Xi %*% ginv(Omegai_hat) %*% t(Xi) %*% M0/T     # T x T
    LM0[i,1] = uhat[i,] %*% PXi %*% uhat[i,]      # scalar

    # Bias correction.
    Hihat = M0 %*% PXi %*% M0     # T x T
    Hihat_diag = diag(Hihat)    # T x 1
    B0[i,1] = (uhat[i,]*uhat[i,]) %*% Hihat_diag    # scalar

    # Variance correction.
    bi_hat = M0 %*% Xi %*% ginv(expm::sqrtm(Omegai_hat))   # T x p
    b_epsilon_i = bi_hat * kronecker(uhat[i,], array(1, c(1,p)))  # T x p
    b_epsilon_i_cumsum = colCumsums(b_epsilon_i)    # T x p
    temp = b_epsilon_i[2:T, ] * b_epsilon_i_cumsum[1:(T-1),]    # (T-1) x p
    if (p==1) {
      V0[i,1] = sum(temp**2)    # scalar
    } else {
      V0[i,1] = sum(rowSums(temp)*rowSums(temp))    # scalar
    }
  }

  LMNT = sum(LM0)
  BNT  = sum(B0)/sqrt(N)
  VNT  = 4*sum(V0)/(N*(T**2))

  # Standardised test statistic
  Jhat = (N**(-1/2)*LMNT-BNT)/sqrt(VNT)

  return(Jhat)
}
