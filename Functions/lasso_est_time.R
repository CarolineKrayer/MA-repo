######################################################################
# Implement C-Lasso estimator for models with time-fixed effects.
# Closely follows implementation in classo package for models without 
# time-fixed effects.
######################################################################
library(Matrix)
library(CVXR)
library(MASS)
library(SparseM)
library(classo)
library(pracma)

source("./Functions/estimators.R")
source("./Functions/demean_data.R")

######## Helper functions from C-Lasso package ###########
# Generate known part of penalty term.
pen_generate = function (b, a, N, p, K, kk) {
  a.out.exp = aperm(array(a, c(K, p, N)), c(3, 2, 1))
  p.temp = (b - a.out.exp)^2
  p.norm = sqrt(apply(p.temp, c(1, 3), sum))
  
  ind = setdiff(1:K, kk)
  pen = apply(as.matrix(p.norm[, ind]), 1, prod)
  
  return(pen)     # N x 1
}

# Convergence criterion.
criterion_time = function (a.old, a.new, b.old, b.new, tol) {
  d = FALSE
  
  a.cri = sum(abs(a.old - a.new))/(sum(abs(a.old)) + 1e-04)
  b.cri = mean(abs(b.old - b.new))/(mean(abs(b.old)) + 1e-04)
  
  if (a.cri < tol & b.cri < tol) {
    d = TRUE
  }
  
  return(d)
}

########################## Mosek solver #################################
# Input:
#       Datamatrix: matrix with regressors stacked for objective function (NT x Np)
#       Y_demean_vector: demeaned dependent variable (NT x 1)
#       penalty: known part of penalty term (N x 1)
#       N: cross-sectional dimension
#       T: time series dimension
#       K: number of groups considered
#       p: number of regressors
#       lambda: tuning parameter for estimation
# Output:
#       beta: individual-specific coefficient estimates (N x p)
#       alpha: group-specific coefficient estimates (K x p)

opt_mosek_time = function (Datamatrix, Y_demean_vector, penalty, N, T, K, p, lambda) {
    # Prerequisites - tolerance criterion  
    prob = list(sense = "min")
    prob$dparam = list(INTPNT_CO_TOL_REL_GAP=1e-5)

    prob$c = c(rep(0, N * (2 * p + T + 2) + p), rep(1/(N * T), N),
        penalty * c(lambda/N))

    # linear constraint: matrix A
    A.y = cbind(Datamatrix, Diagonal(T * N), Matrix(0, T * N,
        N * (p + 4) + p))

    A.0 = cbind(Diagonal(N * p), Matrix(0, N * p, T * N), -Diagonal(N *
        p), -matrix(diag(p), N * p, p, byrow = TRUE), Matrix(0, N * p,
        N * 4))

    A.nhalf = cbind(Matrix(0, N, N * (2 * p + T) + p), Diagonal(N),
        Matrix(0, N, N), -Diagonal(N)/2, Matrix(0, N, N))

    A.phalf = cbind(Matrix(0, N, N * (2 * p + T) + p), Matrix(0, N,
        N), Diagonal(N), -Diagonal(N)/2, Matrix(0, N, N))

    A = rbind(A.y, A.0, A.nhalf, A.phalf)
    prob$A = as(A, "CsparseMatrix")

    # linear constraint: upper bound and lower bound
    prob$bc = rbind(blc = c(Y_demean_vector, rep(0, N * p), rep(-1/2, N), rep(1/2,
        N)), buc = c(Y_demean_vector, rep(0, N * p), rep(-1/2, N), rep(1/2, N)))
    # t_i \geq 0
    prob$bx = rbind(blx = c(rep(-Inf, N * (2 * p + T + 2) + p), rep(0,
        N), rep(-Inf, N)), bux = c(rep(Inf, N * (2 * p + T + 4) + p)))

    # Conic constraints
    CC = list()
    bench = N * (2 * p + T) + p

    for (i in 1:N) {
        s.i = bench + i
        r.i = bench + N + i
        nu.i = (N * p + (i - 1) * T + 1):(N * p + i * T)
        w.i = bench + 3 * N + i
        mu.i = (N * (T + p) + (i - 1) * p + 1):(N * (T + p) + i *
            p)
        CC = cbind(CC, list("QUAD", c(r.i, nu.i, s.i)), list("QUAD",
            c(w.i, mu.i)))
    }
    prob$cones = CC
    rownames(prob$cones) = c("type", "sub")

    # Call mosek solver
    mosek.out = mosek(prob, opts = list(verbose = 0))

    est = mosek.out$sol$itr$xx
    beta = est[1:(N*p)]
    alpha = est[(N*(2*p + T) + 1):(N*(2*p + T) + p)]

    return(list(beta = beta, alpha = alpha))
}

####### Post-Lasso estimation and bias correction. #######
# Input:
#       group.est: array with estimated group number for each individual
#       a.out: matrix with estimated group coefficients (K x p)
#       X: original regressor matrices (N x T x p)
#       Y: matrix with dependent variable (N x T)
#       K: number of groups assumed
#       p: number of regressors
#       N: cross-section dimension
#       T: time-series dimension
#       bias_corr: dummy indicating whether bias correction should be performed
# Output:
#       group_est_corr: (bias corrected) post-Lasso estimates of group-specific 
#                       parameters (K x p)

post_corr_time = function(group.est, a.out, X, Y, K, p, N, T, bias_corr) {
  group_est_corr = matrix(0, K, p)

  for (k in 1:K) {
    # Gather group membership.
    group_k = (group.est == k)
    N_k = sum(group_k)

    if (N_k >= 2 * p/T) {
      indiv = 1:N
      group_index = indiv[group_k]

      # Demean data within group k.
      if (N_k > 1) {
        Y_k = Y[group_index, ]
        X_k = X[group_index, , ]
        
        demean_data = indiv_time_demean(X=X_k, Y=Y_k)
        Y_k_demean_vector = demean_data$Y_demean_vector
        X_k_demean_vector = demean_data$X_demean_vector        
      } else {
        # Special case of one individual only.
        Y_k = array(Y[group_index, ], c(N_k, T))
        X_k = array(X[group_index, , ], c(N_k, T, p))
        
        demean_data = indiv_demean(X=X_k, Y=Y_k)
        Y_k_demean_vector = demean_data$Y_demean_vector
        X_k_demean_vector = demean_data$X_demean_vector  
      }

      # Perform post estimation.
      a_k = lsfit(x=X_k_demean_vector, y=Y_k_demean_vector, intercept=FALSE)$coefficients
        
      # Perform bias correction if specified.
      if (bias_corr == TRUE) {
        bias_k = classo::SPJ_PLS(t=T, y=Y_k_demean_vector, x=X_k_demean_vector)
        group_est_corr[k, ] = 2 * a_k - bias_k        
      } else {
        group_est_corr[k, ] = a_k
      }

    } else {
      group_est_corr[k, ] = a.out[k, ]
    }
  }

  return(group_est_corr)
}


####### C-Lasso estimation for model with time-fixed effects - using CVX. ######
# Input:
#       X: regressor array of dimension N x T x p
#       Y: dependent variable of dimension N x T 
#       X_demean_vector: demeaned regressor in vector format (NT x p)
#       Y_demean_vector: demeaned dependent variable in vector format (NT x 1)
#       Datamatrix: matrix with regressors stacked for objective function (NT x Np)
#       K: number of groups used in estimation
#       lambda: tuning parameter (scalar)
#       beta0: initialised individual-specific coefficients
#       R: maximal number of iterations in optimisation algorithm
#       tol: tolerance level for convergence in optimisation algorithm
#       post_est: dummy indicating whether post estimation is performed
#       bias_corr: dummy indicating whether bias correction is performed
# Output:
#       b.est: estimated individual coefficients (N x p)
#       a.out: estimated group-specific parameters (K x p)
#       group.est: estimated group number for each individual (N x 1)
#       converge: dummy indicating whether convergence criterion is met

lasso_est_time = function(X, Y, X_demean_vector, Y_demean_vector,
                          Datamatrix, K, lambda, beta0=NULL, R=500, tol=1e-04,
                          post_est=TRUE, bias_corr=FALSE) {
  N = dim(X)[1]
  T = dim(X)[2]
  p = dim(X)[3]

  # Initialise coefficient estimates.
  if(is.null(beta0)){
    beta0 = hetero_coef_est(X=X, Y=Y)
  }

  b.out = array(beta0, c(N, p, K))
  a.out = matrix(0, K, p)

  b.old = matrix(1, N, p)
  a.old = matrix(1, 1, p)

  # Optimisation algorithm.
  for (r in 1:R) {
    for (k in 1:K) {
      gamma = pen_generate(b=b.out, a=a.out, N=N, p=p, K=K, kk=k)

      b = Variable(p, N)
      a = Variable(p)
      A = matrix(1, nrow = 1, ncol = N)
      obj1 = t(norm2(b - a %*% A, axis = 2)) %*% gamma

      obj = Minimize(sum_squares(Y_demean_vector - Datamatrix %*% vec(b))/(N * T) +
                       obj1 * (lambda/N))
      Prob = Problem(obj)


      prob_data = get_problem_data(Prob, solver = "ECOS")
      if (packageVersion("CVXR") > "0.99-7") {
        ECOS_dims = ECOS.dims_to_solver_dict(prob_data$data[["dims"]])
      } else {
        ECOS_dims = prob_data$data[["dims"]]
      }
      solver_output = ECOSolveR::ECOS_csolve(c = prob_data$data[["c"]],
                                              G = prob_data$data[["G"]],
                                              h = prob_data$data[["h"]],
                                              dims = ECOS_dims,
                                              A = prob_data$data[["A"]],
                                              b = prob_data$data[["b"]])

      if (packageVersion("CVXR") > "0.99-7") {
        direct_soln = unpack_results(Prob, solver_output, prob_data$chain, prob_data$inverse_data)
      } else {
        direct_soln = unpack_results(Prob, "ECOS", solver_output)
      }

      a.out[k, ] = direct_soln$getValue(a)
      b.out[, , k] = matrix(direct_soln$getValue(b), N, p, byrow = TRUE)
    }

    # Check whether convergence criterion is met.
    a.new = a.out[K, ]
    b.new = b.out[, , K]

    if (criterion_time(a.old=a.old, a.new=a.new, b.old=b.old, b.new=b.new, tol=tol)) {
      break
    }
    # Update coefficient estimates.
    a.old = a.out[K, ]
    b.old = b.out[, , K]
  }

  # Classify individuals to closest group.
  a.out.exp = aperm(array(a.out, c(K, p, N)), c(3, 2, 1))
  d.temp = (b.out - a.out.exp)**2
  dist = sqrt(apply(d.temp, c(1, 3), sum))
  group.est = apply(dist, 1, which.min)

  # Post estimation and bias correction.
  if (post_est) {
    a.out = post_corr_time(group.est = group.est,
                           a.out = a.out,
                           X = X,
                           Y = Y,
                           K = K,
                           p = p,
                           N = N,
                           T = T,
                           bias_corr = bias_corr)
  }

  # Get individual coefficients from group-specific ones.
  b.est = matrix(999, N, p)
  for (i in 1:N) {
    group = group.est[i]
    b.est[i, ] = a.out[group, ]
  }

  return(list(b.est = b.est, 
              a.out = a.out, 
              group.est = group.est,
              converge = (r < R)))
}


## C-Lasso estimation for model with time-fixed effects - using RMosek. ##
# Input:
#       X: regressor array of dimension N x T x p
#       Y: dependent variable of dimension N x T 
#       X_demean_vector: demeaned regressor in vector format (NT x p)
#       Y_demean_vector: demeaned dependent variable in vector format (NT x 1)
#       Datamatrix: matrix with regressors stacked for objective function (NT x Np)
#       K: number of groups used in estimation
#       lambda: tuning parameter (scalar)
#       beta0: initialised individual-specific coefficients
#       R: maximal number of iterations in optimisation algorithm
#       tol: tolerance level for convergence in optimisation algorithm
#       post_est: dummy indicating whether post estimation is performed
#       bias_corr: dummy indicating whether bias correction is performed
# Output:
#       b.est: estimated individual coefficients (N x p)
#       a.out: estimated group-specific parameters (K x p)
#       group.est: estimated group number for each individual (N x 1)
#       converge: dummy indicating whether convergence criterion is met

lasso_est_time_mosek = function(X, Y, X_demean_vector, Y_demean_vector,
                                Datamatrix, K, lambda, beta0=NULL, R=500, tol = 1e-04,
                                post_est=TRUE, bias_corr=FALSE) {
  N = dim(X)[1]
  T = dim(X)[2]
  p = dim(X)[3]

  # Initialise coefficient estimates
  if(is.null(beta0)){
    beta0 = hetero_coef_est(X=X, Y=Y)
  }

  b.out = array(beta0, c(N, p, K))
  a.out = matrix(0, K, p)

  b.old = matrix(1, N, p)
  a.old = matrix(1, 1, p)

  # Optimisation algorithm.
  for (r in 1:R) {
    for (k in 1:K) {
      penalty.out = pen_generate(b=b.out, a=a.out, N=N, p=p, K=K, kk=k)

      mosek.out = opt_mosek_time(Datamatrix=Datamatrix,
                                 Y_demean_vector=Y_demean_vector,
                                 penalty=penalty.out,
                                 N=N,
                                 T=T,
                                 K=K,
                                 p=p,
                                 lambda=lambda)
      a.out[k, ] = mosek.out$alpha
      b.out[, , k] = matrix(mosek.out$beta, N, p, byrow = TRUE)

    }

    # Check whether convergence criterion is met.
    a.new = a.out[K, ]
    b.new = b.out[, , K]

    if (criterion_time(a.old=a.old, a.new=a.new, b.old=b.old, b.new=b.new, tol=tol)) {
      break
    }
    
    # Update coefficient estimates.
    a.old = a.out[K, ]
    b.old = b.out[, , K]
  }

  # Classify individuals to closest group.
  a.out.exp = aperm(array(a.out, c(K, p, N)), c(3, 2, 1))
  d.temp = (b.out - a.out.exp)**2
  dist = sqrt(apply(d.temp, c(1, 3), sum))
  group.est = apply(dist, 1, which.min)


  # Post estimation and bias correction.
  if (post_est == TRUE) {
    a.out = post_corr_time(group.est = group.est,
                           a.out = a.out,
                           X = X,
                           Y = Y,
                           K = K,
                           p = p,
                           N = N,
                           T = T,
                           bias_corr = bias_corr)
  }

  # Get individual coefficients from group-specific ones.
  b.est = matrix(999, N, p)
  for (i in 1:N) {
    group = group.est[i]
    b.est[i, ] = a.out[group, ]
  }

  return(list(b.est = b.est, 
              a.out = a.out, 
              group.est = group.est,
              converge = (r < R)))
}
