################################################################################
#Script-file: zlm_test1.R
#Description: check conditions which make hurdle regression fail to converge
#...fraction of genes active / average transformed count: et
################################################################################

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
  
  #Genes not yet filtered
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
#Try one gene
#Check if it gave an error
test1 <- zlm(
  ~age.group*days_post + (1|indiv),
  sca = sca_l[1,],
  method = "glmer", ebayes = FALSE
)
summary(test1)

################################################################################
#Using zlm SingleCellAssay glmer
#Estimate the proportion of genes which give errors:

#Try for a subset of the genes
n_check <- 500 #dim(sca_l)[1]
pass <- rep(TRUE, n_check)
for (i in 1:n_check) {
  sca_l_1gene <- sca_l[i, ]
  temp <- zlm(
    ~age.group*days_post + (1|indiv),
    sca = sca_l_1gene,
    method = "glmer", ebayes = FALSE
  )

  pass[i] <- (
    as.vector(!is.na(temp@coefC)[[1]]) &
    as.vector(!is.na(temp@coefC)[[2]]) &
    as.vector(!is.na(temp@coefC)[[3]]) &
    as.vector(!is.na(temp@coefC)[[4]]) &
    as.vector(!is.na(temp@coefD)[[1]]) &
    as.vector(!is.na(temp@coefD)[[2]]) &
    as.vector(!is.na(temp@coefD)[[3]]) &
    as.vector(!is.na(temp@coefD)[[4]]) 
  )
}

raw2_gene %<>% mutate("pass" = NA) 
raw2_gene$pass[1:n_check] = pass

summary(raw2_gene$pass)

#####
ggplot(data = raw2_gene %>% dplyr::filter(pass==TRUE), aes(x = frac_active)) +
  geom_histogram(binwidth = 0.05, boundary = 0) +
  coord_cartesian(xlim = c(0,1)) +
  ggtitle("frac_active when zlm passes") +
  theme_bw()

ggplot(data = raw2_gene %>% dplyr::filter(pass==FALSE), aes(x = frac_active)) +
  geom_histogram(binwidth = 0.05, boundary = 0) +
  coord_cartesian(xlim = c(0,1)) +
  ggtitle("frac_active when zlm fails") +
  theme_bw()

quantile(raw2_gene$avg_et, 0.975)

#####
ggplot(data = raw2_gene %>% dplyr::filter(pass==TRUE), aes(x = avg_et)) +
  geom_histogram(binwidth = 0.2, boundary = 0) +
  coord_cartesian(xlim = c(0,quantile(raw2_gene$avg_et, 0.975))) +
  ggtitle("avg_et when zlm passes (upto 97.5%ile of overall avg_et)") +
  theme_bw()

ggplot(data = raw2_gene %>% dplyr::filter(pass==FALSE), aes(x = avg_et)) +
  geom_histogram(binwidth = 0.2, boundary = 0) +
  coord_cartesian(xlim = c(0,quantile(raw2_gene$avg_et, 0.975))) +
  ggtitle("avg_et when zlm fails") +
  theme_bw()

#####
ggplot(data = raw2_gene %>% dplyr::filter(pass==TRUE), aes(x = frac_active, y = avg_et)) +
  stat_bin2d(binwidth = c(0.05, 0.2)) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,quantile(raw2_gene$avg_et, 1))) +
  theme_bw() +
  scale_x_continuous(minor_breaks = seq(0, 1, 0.05), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(minor_breaks = seq(-2, 15, 0.2), breaks = seq(-2, 15, 5)) +
  ggtitle("frac_active, avg_et when zlm passes")

ggplot(data = raw2_gene %>% dplyr::filter(pass==FALSE), aes(x = frac_active, y = avg_et)) +
  stat_bin2d(binwidth = c(0.05, 0.2)) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,quantile(raw2_gene$avg_et, 1))) +
  theme_bw() +
  scale_x_continuous(minor_breaks = seq(0, 1, 0.05), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(minor_breaks = seq(-2, 15, 0.2), breaks = seq(-2, 15, 5)) +
  ggtitle("frac_active, avg_et when zlm fails")

#####
with(
  raw2_gene %>% dplyr::filter(pass==TRUE),
  summary(frac_active)
)

with(
  raw2_gene %>% dplyr::filter(pass==FALSE),
  (summary(frac_active))
)

with(
  raw2_gene %>% dplyr::filter(pass==TRUE),
  summary(avg_et)
)

with(
  raw2_gene %>% dplyr::filter(pass==FALSE),
  (summary(avg_et))
)
