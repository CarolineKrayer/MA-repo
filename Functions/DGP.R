##########################################################
# Implement data generating processes of simulation study.
##########################################################
library(pracma)

# Input:
#       DGP: index of DGP to be considered (integer between 1 and 8 or "power_DGP")
#       time_effect: dummy indicating inclusion of time-fixed effects
#       N: cross section dimension
#       T: time dimension
#       disturbance: disturbance term added to true coefficient
# Output:
#       X: regressor matrices (N x T x p)
#       Y: dependent variable (N x T)
#       beta: true coefficient matrix (N x p)

generate_data = function(DGP, time_effect, N, T, disturbance=NULL){
  # Fix number of regressors.
  p = 2

  # Generate coefficients.
  if (DGP == 1 || DGP == 2 || DGP == 3 || DGP == 4 || DGP == 6) {
    if (DGP == 1 || DGP == 2 || DGP == 3 || DGP == 4){
      N1 = floor(N*0.3)
      N2 = floor(N*0.3)      
    } else {
      N1 = floor(N*0.7)
      N2 = floor(N*0.05)
    }
    beta = array(0, c(N, p))
    beta[1:N1, ] = repmat(c(0.5, -0.5), N1, 1)
    beta[(N1+1):(N1+N2), ] = repmat(c(-0.5, 0.5), N2, 1)
  }

  if (DGP == 5) {
    N1 = floor(N*0.3)
    N2 = floor(N*0.3)      
    
    beta = zeros(N, p)
    beta[1:N1, ] = repmat(c(0.5, -0.5), N1, 1)
    beta[(N1+1):(N1+N2), ] = repmat(c(0.3, -0.2), N2, 1)
  }
  
  if (DGP == 7) {
    beta = zeros(N, p)
    beta[, 1] = rnorm(N, 0.5, 1)
    beta[, 2] = runif(N, min=-0.5, max=0.5)
  }
  
  if (DGP == 8) {
    beta = repmat(c(0.5, -0.5), N, 1)
  }

  if (DGP == "power_DGP") {
    beta = repmat(c(0.5, -0.1), N, 1)
  }
  
  # Add disturbances to coefficients of DGP 4, 8 and "power_DGP".
  if (DGP == 4 || DGP == 8 || DGP == "power_DGP") {
    if (is.null(disturbance)) {
      dstb = 0.1
    } else {
      dstb = disturbance
    }

    beta = beta + dstb * rnorm(N*p, 0, 1)
    
    # Make sure that there is no unit root or non-stationary process.
    for (i in 1:N) {
      if (beta[i, 2]>=1) {
        beta[i, 2] = runif(1, min= 0, max=0.9)
      } else if (beta[i, 2]<=-1) {
        beta[i, 2] = runif(1, min=-0.9, max=0)
      } 
    }
  }
  
  # Generate regressors and dependent variable.
  if (DGP == 1) {
    # Static panel.
    mu = array(rnorm(N, 0, 1), c(N, 1))
    X = array(0, c(N, T, p))
    for (ppp in 1:p) {
      X[, , ppp] = rnorm(N*T, 0, 1) + repmat(mu, 1, T)
    }
    u = array(rnorm(N*T, 0, 1), c(N, T))
    Y = repmat(mu, 1, T) + u
    for (ppp in 1:p) {
      Y = Y + repmat(array(beta[ , ppp], c(N,1)), 1, T) * X[ , , ppp]
    }
    
  } else {
    
    # Dynamic panel.
    TTT = 1000 + T   # Add further time periods which are cancelled later.
    mu = array(rnorm(N, 0, 1), c(N, 1))

    # Generate fixed effects and regressors.
    if (time_effect == TRUE) {
      gamma = array(rnorm(TTT, 0, 1), c(1, TTT))
      X = array(0, c(N, TTT, p))
      X[ , , 1] = rnorm(N*TTT) + repmat(mu, 1, TTT) + repmat(gamma, N, 1)
    } else {
      X = array(0, c(N, TTT, p))
      X[ , , 1] = rnorm(N*TTT) + repmat(mu, 1, TTT)
    }

    # Coefficient of fixed effects.
    coef_indiv_effect = 1
    coef_time_effect = 1
    if (time_effect == TRUE) {
      coef_indiv_effect = 0.5
      coef_time_effect = 0.5
    }

    # Generate error term and dependent variable.
    u = array(rnorm(N*TTT, 0, 1), c(N, TTT))
    Y = zeros(N, TTT)
    Y[ , 1] = zeros(N, 1)
    for (ttt in 2:TTT) {
      X[ , ttt, 2] = Y[, ttt-1]
      Y[ , ttt] = coef_indiv_effect * mu + u[ , ttt]

      if (time_effect == TRUE) {
        Y[ , ttt] = Y[ , ttt] + coef_time_effect * repmat(gamma[ , ttt], N, 1)
      }
      
      for (ppp in 1:p) {
        Y[ , ttt] = Y[ , ttt] + beta[ , ppp] * X[ , ttt, ppp]
      }
    }

    X = X[ , (TTT-T+1):TTT, ]
    Y = Y[ , (TTT-T+1):TTT]
    
  }

  return(list(X=X, Y=Y, beta=beta))
}
