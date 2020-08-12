############################
# Create powerfunction plot.
############################
library(ggplot2)
library(plyr)
library(reshape2)

rm(list=ls())
out_path = "./Results/"

# Set parameters.
N_vec = c(40, 80, 40, 20)
T_vec = c(10, 10, 40, 40)
M = 2
Kmax = 5
alpha = 0.05/Kmax
time_effect = TRUE

# Load data.
data = mapply(function(N, T) readRDS(df, 
                              file=paste0(out_path, "Powerfct_homo_dist_N", N, "_T", T, "_Mrep", M, "_alpha", alpha, "_time_effect", time_effect, ".rds")),
              N_vec,
              T_vec,
              SIMPLIFY=FALSE)

#################################
#### AB HIER GLEICH WEITER!!!
# vllt klappts mit SIMPLIFY=TRUE?!?
#################################
# Merge results in dataframe.
df = merge(data[[1]], data[[2]], by="disturbances")
df = merge(df, data[[3]], by="disturbances")
df = merge(df, data[[4]], by="disturbances")
df_plot = melt(df, id.vars=c("disturbances"))

# Plot the powerfunction.
final = ggplot(df_plot, aes(x=disturbances, y=value, group=variable)) + 
  geom_point(aes(col=variable), alpha=3, size=2) +
  scale_color_manual(values=c("#f8766c", "#54f75f", "#00bec4", "#f7b96e"),
                     name="Sample size",
                     labels=c("N=40, T=10", "N=80, T=10", "N=40, T=40", "N=20, T=40")) +
  labs(
    x = "Standard deviation of disturbance term",
    y = (expression(beta))) + 
  xlim(0, 0.35) +
  geom_hline(yintercept = 0.01, lty=2) +
  geom_text(aes(x=0.3, y=0.04, label='alpha==0.01'), 
            parse=TRUE, size=6) +
  theme_minimal() +
  theme(plot.margin=margin(0.4, 0.4, 0.4, 0.4, "cm"),
        text = element_text(size=20)) +
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1, "cm"),
        legend.position = c(0.85, 0.5),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  theme(legend.title.align = 0.5) + 
  guides(colour = guide_legend(override.aes = list(size=3)))

# Save plot.
ggsave(final, file=paste0(out_path, "Powerfct_homo_dist.png"),
       height=7, width=8)
