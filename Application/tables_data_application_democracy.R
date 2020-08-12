#####################################
# Create tables for data application.
#####################################
library(xtable)
library(rlist)

rm(list=ls())
out_path = "./Application/Results/"

application = "democracy"
Kmax = 5
ccons = c(0.25, 0.5, 0.75)
all_countries = TRUE
n_groups = 2  # Number of groups determined by test.
lambda_pos = 1  # Tuning parameter selected.

# Load estimation results.
results = list.load(file=paste0(out_path, "Results_application", application, "_Kmax", Kmax, "_all_countries", all_countries, ".rdata"))

# Construct table with test statistics and corresponding p-values.
table1 = matrix(NA, nrow = 5, ncol = 6, dimnames = list(c("K=1", "K=2", "K=3", "K=4", "K=5"), c("Test Statistic, c=0.25", "p-Value, c=0.25", "Test statistic, c=0.5", "p-Value, c=0.5", "Test statistic, c=0.75", "p-Value, c=0.75")))
table1[ , c(1, 3, 5)] = results$test_stat
table1[ , c(2, 4, 6)] = results$p_test_stat

table1_tex = xtable(table1, digits=4)
caption(table1_tex) = "Test Statistics with p-Values \n"

print.xtable(x=table1_tex,
             type="latex",
             file=paste0(out_path, "test_stat_table_application_", application, "_all_countries", all_countries, ".tex"),
             caption.placement = "bottom",
             booktabs=TRUE
             )

# Create table with results on information criterion.
table2 = matrix(NA, nrow=6, ncol=3, dimnames=list(c("K=1", "K=2", "K=3", "K=4", "K=5", "K=N"), c("c=0.25", "c=0.5", "c=0.75")))
table2[ , 1] = results$IC[ , 1]
table2[ , 2] = results$IC[ , 2]
table2[ , 3] = results$IC[ , 3]

table2_tex = xtable(table2, digits=3)
caption(table2_tex) = "Values of the Information Criterion \n"
                             
print.xtable(x=table2_tex,
             type="latex",
             file=paste0(out_path, "IC_table_application_", application, "_all_countries", all_countries, ".tex"),
             caption.placement = "bottom",
             booktabs = TRUE
)

# Create table with estimator results.
table3 = matrix(NA, nrow=5, ncol=3*(n_groups+1), dimnames=list(c("Democracy", "Lag 1", "Lag 2", "Lag 3", "Lag 4"), c("Estimate", "SE", "p-Value", "Estimate Group 1", "SE Group 1", "p-Value Group 1", "Estimate Group 2", "SE Group 2", "p-Value Group 2")))

table3[1:5, 1] = results$common_coef_est_corr[1, ]
table3[1:5, 2] = results$se[1, , 1, 1]
table3[1:5, 3] = results$p_estim[1, , 1, 1]

table3[1:5, c(4, 7)] = t(results$post_lasso_est_corr[1:n_groups, , n_groups, lambda_pos])
table3[1:5, c(5, 8)] = t(results$se[1:n_groups, , n_groups, lambda_pos])
table3[1:5, c(6, 9)] = t(results$p_estim[1:n_groups, , n_groups, lambda_pos])

table3_tex = xtable(table3, digits=3)
caption(table3_tex) = "Estimation Results \n"

print.xtable(x=table3_tex,
             type="latex",
             file=paste0(out_path, "estim_table_application_", application, "_all_countries", all_countries, ".tex"),
             caption.placement = "bottom",
             booktabs = TRUE
)



