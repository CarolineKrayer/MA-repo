########################################################################
# Compute empirical rejection frequencies and frequencies of number of 
# groups determined by the test and information criterion, respectively.
########################################################################
library(classo)
library(Rmosek)
library(pbapply)
library(parallel)
library(rlist)
library(xtable)

rm(list=ls())
out_path = "./Simstudy/Results/"

# Compute frequencies.
# Input:
#       DGP: DGP considered (integer between 1 and 8)
#       N: cross-section dimension
#       T: time-series dimension
#       Kmax: maximal number of groups considered
#       M: number of Monte Carlo repetitions
#       alpha: significance level for hypothesis test
#       time_effect: dummy indicating inclusion of time-fixed effects
#       bias_corr: dummy indicating whether bias correction should be performed
#       ccons: constants considered in tuning parameter lambda
# Output:
#       rej_freq: rejection frequencies of hypothesis test (Kmax x length(ccons))
#       est_group_number_freq: frequency of estimated number of groups with test 
#                              (Kmax x length(ccons))
#       info_criterion_freq: frequency of estimated number of groups with 
#                            information criterion (... x length(ccons))

frequencies = function(DGP, N, T, Kmax, M, alpha, time_effect, bias_corr, ccons=c(0.25, 0.5, 0.75)){
  results = list.load(file=paste0(out_path, "Test_IC_results_DGP", DGP, "_time_effect", time_effect, "_N", N, "_T", T, "_Kmax", Kmax, "_Mrep", M, "_bias_corr", bias_corr,".rdata"))
  test_stats = results[1, ]
  info_criterions = results[2, ]

  # Compute empirical rejection frequencies of hypothesis test.
  rej = lapply(test_stats, function(x) as.numeric(x >= qnorm(1-alpha, 0, 1)))
  rej_array = array(as.numeric(unlist(rej)), dim=c(Kmax, length(ccons), M))
  rej_freq = rowSums(rej_array, dims=2)/M

  # Determine estimated number of groups with test.
  est_group_number = apply(rej_array, c(2, 3), function(x) min(which(x == 0)))  # length(ccons) x M
  est_group_number_freq = apply(est_group_number, 1, function(x) table(factor(x, levels=c(1:Kmax, Inf)))/M)    # Kmax x length(ccons)

  # Determine number of groups selected with information criterion.
  info_criterion_array = array(as.numeric(unlist(info_criterions)), dim=c(length(ccons), M)) # length(ccons) x M
  info_criterion_freq = apply(info_criterion_array, 1, function(x) table(factor(x, levels=1:(Kmax+1)))/M)   # Kmax x length(ccons)
  
  return(list(rej_freq=rej_freq,
              est_group_number_freq=est_group_number_freq,
              info_criterion_freq=info_criterion_freq)
  )
}

################################################################################
# Simulation parameters.
N_vec = c(40, 80, 40, 20)
T_vec = c(10, 10, 40, 40)
M = 1000
Kmax = 5
DGP = 6
bias_corr = TRUE
time_effect = TRUE
alpha1 = 0.05
alpha3 = alpha1/Kmax
lambda_pos = 2

# Create LaTex tables with results.
table1 = matrix(NA, nrow = 12, ncol = 9, dimnames = list(c("K=1", "K=2", "K=3", "K=1", "K=2", "K=3", "K=1", "K=2", "K=3", "K=1", "K=2", "K=3"), c("alpha=0.05, c=0.25", "alpha=0.05, c=0.5", "alpha=0.05, c=0.75", "alpha=1/T, c=0.25", "alpha=1/T, c=0.5", "alpha=1/T, c=0.75", "alpha=0.05/Kmax, c=0.25", "alpha=0.05/Kmax, c=0.5", "alpha=0.05/Kmax, c=0.75")))
table2 = matrix(NA, nrow = 12, ncol = 6, dimnames = list(c("N=40, T=10, alpha=0.05", "N=40, T=10, alpha=1/T", "N=40, T=10, alpha=0.05/Kmax", "N=80, T=10, alpha=0.05", "N=80, T=10, alpha=1/T", "N=80, T=10, alpha=0.05/Kmax", "N=40, T=40, alpha=0.05", "N=40, T=40, alpha=1/T", "N=40, T=40, alpha=0.05/Kmax", "N=20, T=40, alpha=0.05", "N=20, T=40, alpha=1/T", "N=20, T=40, alpha=0.05/Kmax"), c("K=1", "K=2", "K=3", "K=4", "K=5", "K=N")))
table3 = matrix(NA, nrow = 4, ncol = 6, dimnames = list(c("N=40, T=10", "N=80, T=10", "N=40, T=40", "N=20, T=40"), c("K=1", "K=2", "K=3", "K=4", "K=5", "K=N")))

for (i in 1:4) {
  alpha2 = 1/T_vec[i]
  
  freq_results1 = frequencies(DGP=DGP, N=N_vec[i], T=T_vec[i], Kmax=Kmax, M=M, alpha=alpha1, time_effect=time_effect, bias_corr=bias_corr)
  freq_results2 = frequencies(DGP=DGP, N=N_vec[i], T=T_vec[i], Kmax=Kmax, M=M, alpha=alpha2, time_effect=time_effect, bias_corr=bias_corr)
  freq_results3 = frequencies(DGP=DGP, N=N_vec[i], T=T_vec[i], Kmax=Kmax, M=M, alpha=alpha3, time_effect=time_effect, bias_corr=bias_corr)

  table1[(1+(i-1)*3):(3+(i-1)*3), c(1:3)] = freq_results1$rej_freq[1:3, ]
  table1[(1+(i-1)*3):(3+(i-1)*3), c(4:6)] = freq_results2$rej_freq[1:3, ]
  table1[(1+(i-1)*3):(3+(i-1)*3), c(7:9)] = freq_results3$rej_freq[1:3, ]
  
  table2[(1+(i-1)*3), ] = freq_results1$est_group_number_freq[1:6, lambda_pos]
  table2[(2+(i-1)*3), ] = freq_results2$est_group_number_freq[1:6, lambda_pos]
  table2[(3+(i-1)*3), ] = freq_results3$est_group_number_freq[1:6, lambda_pos]
  
  table3[(1+(i-1)), ] = t(freq_results1$info_criterion_freq[, lambda_pos])
}

# Table1: Empirical rejection frequencies.
table1_tex = xtable(table1, digits=3)
caption(table1_tex) = paste0("Rejection frequencies DGP", DGP,"time effect", time_effect, "biascorr", bias_corr, "\n")

print.xtable(x=table1_tex,
             type="latex",
             file=paste0(out_path, "rejection_frequencies_DGP", DGP, "_time_effect", time_effect, "_biascorr", bias_corr, "_Mrep", M, "_lambda_pos", lambda_pos,".tex"),
             caption.placement = "bottom",
             booktabs=TRUE
)

# Table 2: Estimated number of groups with testing procedure.
table2_tex = xtable(table2, digits=3)
caption(table2_tex) = paste0("Estimated number of groups with testing procedure DGP", DGP, "time effect", time_effect,  "biascorr", bias_corr, "\n")

print.xtable(x=table2_tex,
             type="latex",
             file=paste0(out_path, "est_group_number_test_DGP", DGP, "_time_effect", time_effect, "_biascorr", bias_corr, "_Mrep", M, "_lambda_pos", lambda_pos,".tex"),
             caption.placement = "bottom",
             booktabs=TRUE
)

# Table 3: Estimated number of groups with information criterion.
table3_tex = xtable(table3, digits=3)
caption(table3_tex) = paste0("Estimated number of groups with information criterion DGP", DGP, "time effect", time_effect, "biascorr", bias_corr, "\n")

print.xtable(x=table3_tex,
             type="latex",
             file=paste0(out_path, "est_group_number_IC_DGP", DGP, "_time_effect", time_effect, "_biascorr", bias_corr, "_Mrep", M, "_lambda_pos", lambda_pos,".tex"),
             caption.placement = "bottom",
             booktabs=TRUE
)