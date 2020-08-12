###########################################
# Further analysis of the estimated groups.
###########################################
library(xtable)
library(rlist)
library(plm)
library(pracma)
library(dplyr)
library(devtools)
library(classo)
library(Rmosek)
library(MASS)
library(matrixStats)
library(Matrix)
library(tidyverse)
library(readstata13)
library(tis)

rm(list=ls())
source("./Functions/estimators.R")
source("./Functions/demean_data.R")
out_path = "./Application/Results/"

# Set parameters.
application = "democracy"
Kmax = 5
ccons = c(0.25, 0.5, 0.75)
all_countries = TRUE
K = 2       # Number of groups determined by test.
lambda_pos = 1     # Tuning parameter selected.

# Merge data as before.
application = "democracy"
data_orig = read.dta13("./Application/DDCGdata_final.dta")
df_orig = data.frame(data_orig)

data = read.dta13("./Application/democracy-balanced-l4.dta")
df = data.frame(data)

df_merge = merge(df, df_orig, by.x=c("id", "year"), by.y=c("wbcode2", "year"), all.x=TRUE)
N = length(unique(df_merge$id))
T = length(unique(df_merge$year))

pdf_merge = pdata.frame(df_merge, index = c("id", "year"), drop.index = FALSE, row.names = TRUE)

# Select countries.
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

# Load results of data application.
results = list.load(file=paste0(out_path, "Results_application", application, "_Kmax", Kmax, "_all_countries", all_countries,".rdata"))

# Load and reshape covariates.
democracy = as.matrix(final$dem.x)[, 5:T]
gdp = as.matrix(final$lgdp)[, 5:T]
gdp_lag1 = as.matrix(plm::lag(final$lgdp, k=1, shift="time"))[, 5:T]
gdp_lag2 = as.matrix(plm::lag(final$lgdp, k=2, shift="time"))[, 5:T]
gdp_lag3 = as.matrix(plm::lag(final$lgdp, k=3, shift="time"))[, 5:T]
gdp_lag4 = as.matrix(plm::lag(final$lgdp, k=4, shift="time"))[, 5:T]
country = as.matrix(final$country_name.x)[, 5:T]

unrest = as.matrix(final$unrest)[, 5:T]
share_invest = as.matrix(final$loginvpc)[, 5:T]/100
tfp = as.matrix(final$ltfp)[, 5:T]/100
reforms = as.matrix(final$marketref)[, 5:T]
trade_share = as.matrix(final$ltrade2)[, 5:T]/100
prim_school = as.matrix(final$lprienr)[, 5:T]/100
sec_school = as.matrix(final$lsecenr)[, 5:T]/100
child_mort = as.matrix(final$lmort)[, 5:T]/100

# Find out group membership of each country.
group_membership = results$group[ , K, ]
group1 = which(group_membership[ , lambda_pos] == 1)
group2 = which(group_membership[ , lambda_pos] == 2)

N1 = length(group1)
N2 = length(group2)

country_group1 = country[group1]
country_group2 = country[group2]
group_memb = list(country_group1, country_group2)
names(group_memb) = c("Group1", "Group2")

final_index_orig = as.matrix(final$id)[, 5:T]
group1_orig = as.numeric(final_index_orig[group1])
group2_orig = as.numeric(final_index_orig[group2])
final_group1 = final[final$id %in% group1_orig, ]
final_group2 = final[final$id %in% group2_orig, ]

# Store group membership in a list.
list.save(group_memb, file=paste0(out_path, "Group_membership_countries.rdata"))

# Gather information of countries in the groups.
gdp_group1 = gdp[group1, ]
gdp_group2 = gdp[group2, ]
gdp_lag1_group1 = gdp_lag1[group1, ]
gdp_lag1_group2 = gdp_lag1[group2, ]
gdp_lag2_group1 = gdp_lag2[group1, ]
gdp_lag2_group2 = gdp_lag2[group2, ]
gdp_lag3_group1 = gdp_lag3[group1, ]
gdp_lag3_group2 = gdp_lag3[group2, ]
gdp_lag4_group1 = gdp_lag4[group1, ]
gdp_lag4_group2 = gdp_lag4[group2, ]
democracy_group1 = democracy[group1, ]
democracy_group2 = democracy[group2, ]

share_invest_group1 = share_invest[group1, ]
share_invest_group2 = share_invest[group2, ]
tfp_group1 = tfp[group1, ]
tfp_group2 = tfp[group2, ]
reforms_group1 = reforms[group1, ]
reforms_group2 = reforms[group2, ]
trade_share_group1 = trade_share[group1, ]
trade_share_group2 = trade_share[group2, ]
unrest_group1 = unrest[group1, ]
unrest_group2 = unrest[group2, ]
prim_school_group1 = prim_school[group1, ]
prim_school_group2 = prim_school[group2, ]
sec_school_group1 = sec_school[group1, ]
sec_school_group2 = sec_school[group2, ]
child_mort_group1 = child_mort[group1, ]
child_mort_group2 = child_mort[group2, ]

# Summary statistics of variables and covariates within groups.
# In levels for easier interpretation.
table = matrix(NA, nrow=10, ncol=4, dimnames=list(c("GDP per capita", "Democracy","Share of investment in GDP", "Trade share of GDP", "TFP", "Economic reforms", "Primary school enrollment", "Secondary school enrollment", "Child mortality", "Social unrest dummy"), c("Mean G1", "Sd G1", "Mean G2", "Sd G2")))

covars_group1 = list(exp(gdp_group1), democracy_group1, exp(share_invest_group1), exp(trade_share_group1), exp(tfp_group1), reforms_group1, exp(prim_school_group1), exp(sec_school_group1), exp(child_mort_group1), unrest_group1)
covars_group2 = list(exp(gdp_group2), democracy_group2, exp(share_invest_group2), exp(trade_share_group2), exp(tfp_group2), reforms_group2, exp(prim_school_group2), exp(sec_school_group2), exp(child_mort_group2), unrest_group2)
table[, 1] = sapply(covars_group1, function(x) mean(matrixcalc::vec(x), na.rm=TRUE))
table[, 2] = sapply(covars_group1, function(x) sd(matrixcalc::vec(x), na.rm=TRUE))
table[, 3] = sapply(covars_group2, function(x) mean(matrixcalc::vec(x), na.rm=TRUE))
table[, 4] = sapply(covars_group2, function(x) sd(matrixcalc::vec(x), na.rm=TRUE))

table_tex = xtable(table, digits=2)
caption(table_tex) = "Descriptive statistics within groups. \n"
label(table_tex) = "tab: descriptive_stats"

print.xtable(x=table_tex,
             type="latex",
             file=paste0(out_path, "descriptive_stats_application_", application, "_all_countries", all_countries,".tex"),
             caption.placement = "bottom",
             booktabs = TRUE
)


# Compute correlation of GDP and democracy within the groups.
cor_group1 = cor(matrixcalc::vec(gdp_group1), matrixcalc::vec(democracy_group1))
cor_group2 = cor(matrixcalc::vec(gdp_group2), matrixcalc::vec(democracy_group2))


# Create variables for transitions to democracies and reversals to non-democracies.
trans_democ = matrix(0, nrow=N, ncol=T-4)
trans_democ = pracma::Reshape(mapply(function(x, y) { as.numeric(x-y==1)}, 
                                     democracy, 
                                     cbind(rep(99, N), democracy[, 1:(T-5)])),
                              N, T-4)
trans_democ_group1 = trans_democ[group1, ]  
trans_democ_group2 = trans_democ[group2, ]

rev_democ = matrix(0, nrow=N, ncol=T-4)
rev_democ = pracma::Reshape(mapply(function(x, y) { as.numeric(x-y==-1)}, 
                                   democracy, 
                                   cbind(rep(99, N), democracy[, 1:(T-5)])),
                            N, T-4)
rev_democ_group1 = rev_democ[group1, ]
rev_democ_group2 = rev_democ[group2, ]

# How many countries are always or never democratic?
democracy_years_group1 = rowSums(democracy_group1)
democracy_years_group2 = rowSums(democracy_group2)

always_democ_group1 = sum(democracy_years_group1==19)/N1
always_democ_group2 = sum(democracy_years_group2==19)/N2

never_democ_group1 = sum(democracy_years_group1==0)/N1
never_democ_group2 = sum(democracy_years_group2==0)/N2


# Number of democratisations and reversals to democracy.
num_trans_democ = sum(rowSums(trans_democ))
num_trans_democ_group1 = sum(rowSums(trans_democ_group1))
num_trans_democ_group2 = sum(rowSums(trans_democ_group2))

num_rev_democ = sum(rowSums(rev_democ))
num_rev_democ_group1 = sum(rowSums(rev_democ_group1))
num_rev_democ_group2 = sum(rowSums(rev_democ_group2))

