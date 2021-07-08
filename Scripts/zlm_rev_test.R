################################################################################
#Script-file: zlm_rev_test.R
#Description: check regressions with coding reversed
#...for age (1 young 0 old) 
#...day post vaccine (1 before 0 after)
################################################################################

#Load saved results of setup1a, setup1b
#...with reverse coding included
#(1 young, 0 old)
#(1 before, 0 after)

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

#Genes filtered by frac_active and avg_et
load(file = "./Data/setup1_b3_fltr_rev_sav.RData")

################################################################################
#Use all cores for faster computation?
if (FALSE) {
  core_lim <- detectCores()
  Sys.setenv("MC_CORES" = core_lim)
  library(parallel)
  options(mc.cores = core_lim)
}
print("Environment number of cores: ")
print(Sys.getenv("MC_CORES"))
print("-----")
print("Options number of cores: ")
print(options()[["mc.cores"]])

#Save memory by only keeping necessary variables
rm(gse, raw, raw_gene, egENSEMBL_tb, egSymbol_tb, gene_info)
rm(core_lim, opt_temp)

################################################################################
#Using zlm SingleCellAssay glmer
#Try, say, just the first 200 genes
#Make a vector marking if in the first 200 (primerid order)
raw2_gene_cut <- raw2_gene %>% dplyr::arrange(primerid) %>% dplyr::slice(1:200)
first_g_lgl <- raw2_gene$primerid %in% raw2_gene_cut$primerid
sum(first_g_lgl)

sca_l_rev <- sca_l_rev[first_g_lgl,]
dim(sca_l_rev)

########################################
#Using both vaccine (before/after) and age as covariates
zlm_glmer_rev0 <- zlm(
  ~age.group*days_post + (1|indiv),
  sca = sca_l_rev,
  method = "glmer", ebayes = FALSE, parallel = TRUE
)

zlm_glmer_rev1 <- zlm(
  ~young*before + (1|indiv),
  sca = sca_l_rev,
  method = "glmer", ebayes = FALSE, parallel = TRUE
)

est_rev0 <- summary(zlm_glmer_rev0, doLRT = TRUE)
est_rev1 <- summary(zlm_glmer_rev1, doLRT = TRUE)

########################################
#Using only vaccine (before/after) as a covariate
zlm_glmer_rev0b <- zlm(
  ~days_post + (1|indiv),
  sca = sca_l_rev,
  method = "glmer", ebayes = FALSE, parallel = TRUE
)

zlm_glmer_rev1b <- zlm(
  ~before + (1|indiv),
  sca = sca_l_rev,
  method = "glmer", ebayes = FALSE, parallel = TRUE
)

est_rev0b <- summary(zlm_glmer_rev0b, doLRT = TRUE)
est_rev1b <- summary(zlm_glmer_rev1b, doLRT = TRUE)

################################################################################

if (FALSE) {
  save.image(file = "./Data/zlm_rev_mid_sav.RData")
}

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
  
  load(file = "./Data/zlm_rev_mid_sav.RData") 
  rm(sca_d00, sca_d42, raw2_samp, raw2_mx_ct, raw2_mx_et)
}

################################################################################
#Convenience for later on:
#...define a function which converts a z-score to a 2-sided-p-value
norm_2sp <- function(x) {
  1 - (pnorm(q = abs(x), mean = 0, sd = 1) - pnorm(q = -abs(x), mean = 0, sd = 1))
}

################################################################################
#Vaccine (before and after) and age as covariates
################################################################################
#Standard coding
est_rev0_dt <- as_tibble(est_rev0[["datatable"]]) %>% 
  dplyr::arrange(primerid, contrast, component)
est_rev0_dt$component %<>% as.factor()
est_rev0_dt$contrast %<>% as.factor()

#Reverse coding
est_rev1_dt <- as_tibble(est_rev1[["datatable"]]) %>% 
  dplyr::arrange(primerid, contrast, component)
est_rev1_dt$component %<>% as.factor()
est_rev1_dt$contrast %<>% as.factor()

#Included genes (first 200 in numerical ENSG order)
raw2_gene_cut %<>% dplyr::arrange(primerid)

#From the test result objects, extract the logFC, D, C part of the genes
#make p-value from z-score, make negative log 10 of p-value

################################################################################
################################################################################
#...interaction (young vaccinated to old vaccinated)
#####
#fold change
inter_rev0_f_tb <- as_tibble(est_rev0[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "age.group67-86yo:days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
inter_rev0_f_tb$fdr <- p.adjust(inter_rev0_f_tb$p_val, method = "BH")
inter_rev0_f_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#discrete
inter_rev0_d_tb <- as_tibble(est_rev0[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "D" & contrast == "age.group67-86yo:days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
inter_rev0_d_tb$fdr <- p.adjust(inter_rev0_d_tb$p_val, method = "BH")
inter_rev0_d_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#continuous
inter_rev0_c_tb <- as_tibble(est_rev0[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "C" & contrast == "age.group67-86yo:days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
inter_rev0_c_tb$fdr <- p.adjust(inter_rev0_c_tb$p_val, method = "BH")
inter_rev0_c_tb %<>% relocate("fdr", .after = "neg_log_10_p")

########################################
########################################
#...interaction with reverse coding
#####
#fold change
inter_rev1_f_tb <- as_tibble(est_rev1[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "young1:before1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
inter_rev1_f_tb$fdr <- p.adjust(inter_rev1_f_tb$p_val, method = "BH")
inter_rev1_f_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#discrete
inter_rev1_d_tb <- as_tibble(est_rev1[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "D" & contrast == "young1:before1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
inter_rev1_d_tb$fdr <- p.adjust(inter_rev1_d_tb$p_val, method = "BH")
inter_rev1_d_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#continuous
inter_rev1_c_tb <- as_tibble(est_rev1[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "C" & contrast == "young1:before1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
inter_rev1_c_tb$fdr <- p.adjust(inter_rev1_c_tb$p_val, method = "BH")
inter_rev1_c_tb %<>% relocate("fdr", .after = "neg_log_10_p")

################################################################################
################################################################################
#...age normal coding
#####
#fold change
age_rev0_f_tb <- as_tibble(est_rev0[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "age.group67-86yo") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
age_rev0_f_tb$fdr <- p.adjust(age_rev0_f_tb$p_val, method = "BH")
age_rev0_f_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#discrete
age_rev0_d_tb <- as_tibble(est_rev0[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "D" & contrast == "age.group67-86yo") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
age_rev0_d_tb$fdr <- p.adjust(age_rev0_d_tb$p_val, method = "BH")
age_rev0_d_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#continuous
age_rev0_c_tb <- as_tibble(est_rev0[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "C" & contrast == "age.group67-86yo") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
age_rev0_c_tb$fdr <- p.adjust(age_rev0_c_tb$p_val, method = "BH")
age_rev0_c_tb %<>% relocate("fdr", .after = "neg_log_10_p")

########################################
########################################
#...age reverse coding
#####
#fold change
age_rev1_f_tb <- as_tibble(est_rev1[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "young1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
age_rev1_f_tb$fdr <- p.adjust(age_rev1_f_tb$p_val, method = "BH")
age_rev1_f_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#discrete
age_rev1_d_tb <- as_tibble(est_rev1[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "D" & contrast == "young1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
age_rev1_d_tb$fdr <- p.adjust(age_rev1_d_tb$p_val, method = "BH")
age_rev1_d_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#continuous
age_rev1_c_tb <- as_tibble(est_rev1[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "C" & contrast == "young1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
age_rev1_c_tb$fdr <- p.adjust(age_rev1_c_tb$p_val, method = "BH")
age_rev1_c_tb %<>% relocate("fdr", .after = "neg_log_10_p")

################################################################################
################################################################################
#...vaccine normal coding
#####
#fold change
vacc_rev0_f_tb <- as_tibble(est_rev0[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev0_f_tb$fdr <- p.adjust(vacc_rev0_f_tb$p_val, method = "BH")
vacc_rev0_f_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#discrete
vacc_rev0_d_tb <- as_tibble(est_rev0[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "D" & contrast == "days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev0_d_tb$fdr <- p.adjust(vacc_rev0_d_tb$p_val, method = "BH")
vacc_rev0_d_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#continuous
vacc_rev0_c_tb <- as_tibble(est_rev0[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "C" & contrast == "days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev0_c_tb$fdr <- p.adjust(vacc_rev0_c_tb$p_val, method = "BH")
vacc_rev0_c_tb %<>% relocate("fdr", .after = "neg_log_10_p")

########################################
########################################
#...vaccine reverse coding
#####
#fold change
vacc_rev1_f_tb <- as_tibble(est_rev1[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "before1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev1_f_tb$fdr <- p.adjust(vacc_rev1_f_tb$p_val, method = "BH")
vacc_rev1_f_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#discrete
vacc_rev1_d_tb <- as_tibble(est_rev1[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "D" & contrast == "before1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev1_d_tb$fdr <- p.adjust(vacc_rev1_d_tb$p_val, method = "BH")
vacc_rev1_d_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#continuous
vacc_rev1_c_tb <- as_tibble(est_rev1[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "C" & contrast == "before1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev1_c_tb$fdr <- p.adjust(vacc_rev1_c_tb$p_val, method = "BH")
vacc_rev1_c_tb %<>% relocate("fdr", .after = "neg_log_10_p")

################################################################################
################################################################################
#From after young to after old
#Using the normal coding (x axis)
#vs
#Using the reverse coding (y axis)
#c: continuous, d: discrete, f: log fold change

vacc_age_comp_tb <- tibble(
  "rev_n" = as.vector(age_rev0_f_tb$coef + inter_rev0_f_tb$coef), 
  "rev_y" = as.vector(age_rev1_f_tb$coef)
)

ggplot(data = vacc_age_comp_tb) +
  geom_point(aes(x=rev_n,y=rev_y)) +
  geom_abline(slope=-1, intercept=0, colour="red")

#
#
#
#
################################################################################
#Vaccine (before and after) only as a covariate
################################################################################

################################################################################
#Standard coding
est_rev0b_dt <- as_tibble(est_rev0b[["datatable"]]) %>% 
  dplyr::arrange(primerid, contrast, component)
est_rev0b_dt$component %<>% as.factor()
est_rev0b_dt$contrast %<>% as.factor()

#Reverse coding
est_rev1b_dt <- as_tibble(est_rev1b[["datatable"]]) %>% 
  dplyr::arrange(primerid, contrast, component)
est_rev1b_dt$component %<>% as.factor()
est_rev1b_dt$contrast %<>% as.factor()

#Included genes (first 200 in numerical ENSG order)
raw2_gene_cut %<>% dplyr::arrange(primerid)

#From the test result objects, extract the logFC, D, C part of the genes
#make p-value from z-score, make negative log 10 of p-value

########################################
#...vaccine normal coding
#####
#fold change
vacc_rev0b_f_tb <- as_tibble(est_rev0b[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev0b_f_tb$fdr <- p.adjust(vacc_rev0b_f_tb$p_val, method = "BH")
vacc_rev0b_f_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#discrete
vacc_rev0b_d_tb <- as_tibble(est_rev0b[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "D" & contrast == "days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev0b_d_tb$fdr <- p.adjust(vacc_rev0b_d_tb$p_val, method = "BH")
vacc_rev0b_d_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#continuous
vacc_rev0b_c_tb <- as_tibble(est_rev0b[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "C" & contrast == "days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev0b_c_tb$fdr <- p.adjust(vacc_rev0b_c_tb$p_val, method = "BH")
vacc_rev0b_c_tb %<>% relocate("fdr", .after = "neg_log_10_p")

########################################
#...vaccine reverse coding
#####
#fold change
vacc_rev1b_f_tb <- as_tibble(est_rev1b[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "before1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev1b_f_tb$fdr <- p.adjust(vacc_rev1b_f_tb$p_val, method = "BH")
vacc_rev1b_f_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#discrete
vacc_rev1b_d_tb <- as_tibble(est_rev1b[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "D" & contrast == "before1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev1b_d_tb$fdr <- p.adjust(vacc_rev1b_d_tb$p_val, method = "BH")
vacc_rev1b_d_tb %<>% relocate("fdr", .after = "neg_log_10_p")

#####
#continuous
vacc_rev1b_c_tb <- as_tibble(est_rev1b[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "C" & contrast == "before1") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene_cut[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_rev1b_c_tb$fdr <- p.adjust(vacc_rev1b_c_tb$p_val, method = "BH")
vacc_rev1b_c_tb %<>% relocate("fdr", .after = "neg_log_10_p")

################################################################################
#Vaccine effect comparison
#Using the normal coding (x axis)
#vs
#Using the reverse coding (y axis)
#c: continuous, d: discrete, f: log fold change

#####

ggplot(tibble(
  "rev_n" = as.vector(vacc_rev0b_c_tb$coef), 
  "rev_y" = as.vector(vacc_rev1b_c_tb$coef))
) +
  geom_point(aes(x=rev_n,y=rev_y))

ggplot(tibble(
  "rev_n" = as.vector(vacc_rev0b_c_tb$p_val), 
  "rev_y" = as.vector(vacc_rev1b_c_tb$p_val))
) +
  geom_point(aes(x=rev_n,y=rev_y))

#####

ggplot(tibble(
  "rev_n" = as.vector(vacc_rev0b_d_tb$coef), 
  "rev_y" = as.vector(vacc_rev1b_d_tb$coef))
) +
  geom_point(aes(x=rev_n,y=rev_y))

ggplot(tibble(
  "rev_n" = as.vector(vacc_rev0b_d_tb$p_val), 
  "rev_y" = as.vector(vacc_rev1b_d_tb$p_val))
) +
  geom_point(aes(x=rev_n,y=rev_y))

#####

ggplot(tibble(
  "rev_n" = as.vector(vacc_rev0b_f_tb$coef), 
  "rev_y" = as.vector(vacc_rev1b_f_tb$coef))
) +
  geom_point(aes(x=rev_n,y=rev_y))

ggplot(tibble(
  "rev_n" = as.vector(vacc_rev0b_f_tb$p_val), 
  "rev_y" = as.vector(vacc_rev1b_f_tb$p_val))
) +
  geom_point(aes(x=rev_n,y=rev_y))
