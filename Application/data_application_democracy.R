###########################################################################
# Perform data application of growth and democracy (Acemoglu et. al. 2019)
###########################################################################
library(readstata13)

rm(list=ls())
source("./Application/data_application.R")
out_path = "./Application/Results/"

# Merge original data of Acemoglu et al (2019) with balanced subpanel.
application = "democracy"
data_orig = read.dta13("./Application/Data/DDCGdata_final.dta")
df_orig = data.frame(data_orig)

data = read.dta13("./Application/Data/democracy-balanced-l4.dta")
df = data.frame(data)

df_merge = merge(df, df_orig, by.x=c("id", "year"), by.y=c("wbcode2", "year"), all.x=TRUE)
N = length(unique(df_merge$id))
T = length(unique(df_merge$year))

pdf_merge = pdata.frame(df_merge, index = c("id", "year"), drop.index = FALSE, row.names = TRUE)

# Select countries.
all_countries = FALSE

if (all_countries==TRUE) {
  final = pdf_merge
} else {
  # Drop countries whose measure of democracy remains constant.
  dem = as.matrix(pdf_merge$dem.x)[, 5:T]
  index_orig = as.matrix(pdf_merge$id)[, 5:T]
  const = rowMeans(dem)
  drop_const1 = which(const == 0)
  drop_const2 = which(const == 1)
  drop_const = c(drop_const1, drop_const2)
  indiv = 1:N
  remain_index = indiv[-drop_const]
  N = length(remain_index)
  
  remain_index_orig = as.numeric(index_orig[remain_index])
  final = pdf_merge[pdf_merge$id %in% remain_index_orig, ]
}

# Form regressor matrix. 
democracy = as.matrix(final$dem.x)[, 5:T]
gdp = as.matrix(final$lgdp)[, 5:T]
lag1 = as.matrix(plm::lag(final$lgdp, k=1, shift="time"))[, 5:T]
lag2 = as.matrix(plm::lag(final$lgdp, k=2, shift="time"))[, 5:T]
lag3 = as.matrix(plm::lag(final$lgdp, k=3, shift="time"))[, 5:T]
lag4 = as.matrix(plm::lag(final$lgdp, k=4, shift="time"))[, 5:T]

T = T-4
p = 5
X = array(0, c(N, T, p))
X[, , 1] = democracy
X[, , 2] = lag1
X[, , 3] = lag2
X[, , 4] = lag3
X[, , 5] = lag4
Y = gdp

# Set parameters for C-Lasso.
Kmax = 5
ccons = c(0.25, 0.5, 0.75)
time_effect = TRUE
bias_corr = TRUE

# Call function for parameter estimation and hypothesis test.
output = data_application(X=X, 
                          Y=Y, 
                          N=N, 
                          T=T, 
                          p=p, 
                          Kmax=Kmax, 
                          time_effect=time_effect,
                          bias_corr = bias_corr,
                          ccons=ccons)

# Save results.
list.save(output, file=paste0(out_path, "Results_application", application, "_Kmax", Kmax, "_all_countries", all_countries, ".rdata"))
