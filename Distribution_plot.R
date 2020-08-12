################################################
# Plot density of test statistic under the null.
################################################
library(ggplot2)
library(cowplot)
library(reshape2)
library(rlist)

rm(list=ls())
out_path = "./Results/"

# Plot density.
# Input:
#       N: cross-section dimension
#       T: time series dimension
#       M: number of Monte Carlo repetitions
#       DGP: data generating process considered (integer between 1 and 8)
#       Kmax: maximal number of groups considered
#       K_true: true number of groups
#       time_effect: dummy indicating inclusion of time-fixed effects
#       lambda_pos: tuning parameter taken into account (position in vector ccons)
#       ccons: all possible tuning parameters
# Output:
#       final_plot: density plot for one DGP and sample size

dist_plot = function(N, T, M, DGP, Kmax, K_true, time_effect, lambda_pos, ccons=c(0.25, 0.5, 0.75)) {
  # Take standard normal distribution as comparison.
  set.seed(864)
  std_norm = rnorm(1000, 0, 1)

  # Load results with and without bias correction. 
  results = list.load(file=paste0(out_path, "Test_IC_results_DGP", DGP, "_time_effect", time_effect, "_N", N, "_T", T, "_Kmax", Kmax, "_Mrep", M, "_bias_corrFALSE", ".rdata"))
  results_bias_corr = list.load(file=paste0(out_path, "Test_IC_results_DGP", DGP, "_time_effect", time_effect, "_N", N, "_T", T, "_Kmax", Kmax, "_Mrep", M, "_bias_corrTRUE", ".rdata"))
  
  test_stats = results[1, ]
  test_stats_bias_corr = results_bias_corr[1, ]
  test_stats_null = matrix(ncol=length(ccons), nrow=M)
  test_stats_null_bias_corr = matrix(ncol=length(ccons), nrow=M)
  for(m in 1:M){
    test_stats_null[m, ] = test_stats[[m]][K_true, ]
    test_stats_null_bias_corr[m, ] = test_stats_bias_corr[[m]][K_true, ]
  } 
  
  # Create dataframe for plotting.
  df_test_stats_null = cbind(test_stats_null[ , lambda_pos], 
                             test_stats_null_bias_corr[ , lambda_pos],
                             std_norm)
  df_test_stats_null = data.frame(df_test_stats_null)
  df_plot <- melt(df_test_stats_null)
  
  # Plot densities.
  final_plot = ggplot(df_plot, aes(x=value, group=variable)) +
                geom_density(aes(fill=variable, color=variable), alpha=.7, size=1.05) +
                geom_hline(yintercept=0, color="#ebebeb", size=1.05) +
                scale_fill_manual(values=c("#f8766c", "#54f75f", "#00bec4"),
                                  name="",
                                  labels=c("Without bias correction", "With bias correction", "Limiting distribution")) +
                scale_color_manual(values=c("#f8766c", "#54f75f", "#00bec4")) +
                labs(
                  x = "Test statistic",
                  y = ("Density")) + 
                theme_minimal() +
                theme(plot.margin=margin(0, 0.5, 0.5, 0.5, "cm"),
                      text=element_text(size=28)) + 
                guides(color=FALSE)
  
  return(final_plot)
}
################################################################################
# Simulation parameters.
N_vec = c(rep(40, 4))
T_vec = rep(c(10, 40), times=2)
DGP_vec = rep(c(5, 42), each=2)
time_effect = TRUE
M = 1000
Kmax = 5
K_true = 3
lambda_pos = 2
ccons = c(0.25, 0.5, 0.75)

plots = mapply(function(n, t, dgp) dist_plot(N=n, 
                                             T=t, 
                                             M=M, 
                                             DGP=dgp, 
                                             Kmax=Kmax, 
                                             K_true=K_true, 
                                             time_effect=time_effect,
                                             lambda_pos=lambda_pos, 
                                             ccons=ccons),
               N_vec,
               T_vec,
               DGP_vec,
               SIMPLIFY=FALSE)

final_plot = plot_grid(plots[[1]] + theme(legend.position="none"), 
                       plots[[2]] + theme(legend.position="none"), 
                       plots[[3]] + theme(legend.position="none"), 
                       plots[[4]] + theme(legend.position="none"),
                       align = "vh",
                       labels = c("A", "B", "C", "D"),
                       label_size=30,
                       nrow = 2)

# Add common legend to the plot. 
legend_b = get_legend(
    plots[[1]] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", panel.spacing = unit(1, "cm")) + 
    guides(color=FALSE)
    )

final = plot_grid(final_plot, legend_b, ncol = 1, rel_heights = c(1, .1))

# Save plot.
ggsave(final, file=paste0(out_path, "Density_plots.png"),
       height=10, width=15)
