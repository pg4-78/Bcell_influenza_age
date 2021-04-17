#Clear
#...variables
rm(list=ls())
#...console
cat("\014\n")
#...graphs
tryCatch(dev.off(), error = function(e) {NULL})
dev.new()

library(tidyverse)
library(magrittr)
library(MAST)

#Load data results of "setup1a.R", "setup1b.R"?
#Alternatively, re-run the setup scripts
if (TRUE) {
  load(file = "./Data/setup1_a_sav.RData")
}

#Save memory by only keeping necessary variables
# rm(list=setdiff(ls(), c("sca_l", "sca_d42", "sca_d00")))

################################################################################
#Convenience for later on:
#...define a function which converts a z-score to a 2-sided-p-value
norm_2sp <- function(x) {
  1 - (pnorm(q = abs(x), mean = 0, sd = 1) - pnorm(q = -abs(x), mean = 0, sd = 1))
}

norm_2sp(0)
norm_2sp(1.96)

################################################################################
#A matrix of individual:gene specific baselines in et = log2(TPM+1)
#Subset with day 0 only
samp_d00_mx <- as_tibble(sca_d00@assays@data@listData[["et"]])

#...Check how many individuals there are
length(levels(raw2_samp$indiv))

#The rows are genes: keep the same
#The columns are individuals
#For every individual: make a column which is the average of their et

#...A template of zeros with genes as rows and individuals as columns
indiv_d00_mx <- as_tibble(
  matrix(0, nrow = dim(samp_d00_mx)[1], ncol = length(levels(raw2_samp$indiv))),
  .name_repair = "unique"
)
colnames(indiv_d00_mx) <- levels(raw2_samp$indiv)

#Get the all the sample et columns corresponding to 1 individual
#Average them
#Put this column into the matrix of baseline averages (columns: individuals)
#...using the previously made template of zeros
for (i in levels(raw2_samp$indiv)) {
  temp_indiv <- samp_d00_mx[,sca_d00@colData@listData[["indiv"]]==i]
  temp_ave_indiv <- tibble("a" = base::rowMeans(temp_indiv))
  indiv_d00_mx[i] <- temp_ave_indiv
}

#...A template of zeros with genes as rows 
#baseline et corresponding to samples at day 42 as columns

base_match_d42_mx <- as_tibble(
  matrix(0, nrow = dim(sca_d42)[1], ncol = dim(sca_d42)[2]),
  .name_repair = "unique"
)
colnames(base_match_d42_mx) <- sca_d42@colData@rownames

#Put copies of these columns to match the order of the day 42 et matrix
#For each individual
#...Go through every column in the collection of day 42 samples
#test if the column belongs to the current individual
#if yes, place a copy of their baseline column into the baseline matrix
for (i in sca_d42@colData@listData[["indiv"]]) {
  base_col <- indiv_d00_mx[,i]
  match_cols <- (sca_d42@colData@listData[["indiv"]] == i)
  for (j in 1:dim(sca_d42)[2]) {
    if (match_cols[j] == TRUE) {
      base_match_d42_mx[,j] <- base_col
    }
  }
}



################################################################################
#Model
zlm_a <- zlm(~age.group, sca_d42)
#(will also need indiv average of day0 in et and vaccine status)

#Save point immediately after the zlm model gets fit (to save time)
if (FALSE) {
  save.image(file = "./Data/zlm1_a_sav.RData")
}

if (FALSE) {
  load(file = "./Data/zlm1_a_sav.RData")
}

#Results
est_a <- summary(zlm_a)

#Table part of results
tb_a <- as_tibble(est_a[["datatable"]]) %>% 
  mutate(p_val = norm_2sp(z))

summary(as_factor(tb_a$component))
summary(as_factor(tb_a$contrast))

#Test genes for the age.group coefficient
test_age <- summary(zlm_a, doLRT = "age.group67-86yo")
test_age2 <- as_tibble(test_age[["datatable"]]) 
test_age3 <- dplyr::filter(test_age2, component == "logFC")
test_age4 <- test_age3 %>% mutate(p_val = norm_2sp(z)) %>% 
  mutate()

ggplot(data = test_age4, aes(x=p_val)) +
  geom_histogram(binwidth = 0.05, boundary = 0)

print(test_age, n = 10, by = "D")

################################################################################
#Gene enrichment
