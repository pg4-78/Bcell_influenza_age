#Immediately after setup1a, setup1b
#Alternatively, load saved results of setup1a, setup1b
if (TRUE) {
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
  library(MAST) #renv::install("RGLab/MAST")
  
  load(file = "./Data/setup1_b3_sav.RData")
}

#Save memory by only keeping necessary variables
# rm(list=setdiff(ls(), c("sca_l", "sca_d42", "sca_d00")))
rm(gse, raw, raw_gene, egENSEMBL_tb, egSymbol_tb, gene_info)

################################################################################
#Checks
if (FALSE) {
  sum(sca_d42@colData@listData[["treatment"]]=="vaccine")
  sum(sca_d42@colData@listData[["treatment"]]=="control")
  
  sum(sca_d00@colData@listData[["treatment"]]=="vaccine")
  sum(sca_d00@colData@listData[["treatment"]]=="control")
  
  #####
  summary(sca_d00@colData@listData[["indiv"]])
  summary(sca_d42@colData@listData[["indiv"]])
  
  summary(sca_l@colData@listData[["indiv"]])
}

################################################################################
#Convenience for later on:
#...define a function which converts a z-score to a 2-sided-p-value
norm_2sp <- function(x) {
  1 - (pnorm(q = abs(x), mean = 0, sd = 1) - pnorm(q = -abs(x), mean = 0, sd = 1))
}

norm_2sp(0)
norm_2sp(1.96)

################################################################################
#Method 1: using baseline et matrix: columns match day42 samples
#...but entries have average day0 et of the corresponding individual

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

####################
#Cycle through each gene
#Perform zlm
#Store the results in a list
#Before attempting loop, first just try one gene
i <-  1
tb_42_v_00 <- tibble(
  et42 = as.matrix(sca_d42@assays@data@listData[["et"]][i,]), 
  age.group = sca_d42@colData@listData[["age.group"]], 
  days_post = sca_d42@colData@listData[["days_post"]],
  CDR = sca_d42@colData@listData[["CDR"]],
  et00 = t(base_match_d42_mx[i,])
)

zlm_out_a <- zlm(
  et42 ~ age.group + et00 + CDR, 
  sca = tb_42_v_00
)

lm_out_a <- lm(
  et42 ~ age.group + et00 + CDR, 
  data = tb_42_v_00
)

################################################################################
#Gene enrichment
