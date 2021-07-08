################################################################################
#Script-file: zlm1a.R
#Description: Main hurdle model for vaccine/age part B
#... volcano plots, gene enrichment, p-value distribution histograms
################################################################################

#Immediately after zlm1a
#or alternatively, load save
#or go to the save after the summary zlm command
if (FALSE) {
  #Clear
  #...variables
  rm(list=ls())
  #...console
  cat("\014\n")
  #...graphs
  tryCatch(dev.off(), error = function(e) {NULL})
  dev.new()
  
  library(tidyverse)
  library(ggrepel) #for volcano plot gene labels
  library(magrittr)
  library(MAST) #renv::install("RGLab/MAST")
  
  #Genes filtered by frac_active and avg_et
  load(file = "./Data/zlm1a_sav.RData")
}

################################################################################
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

################################################################################
#Check
summary(zlm_glmer@converged)
 
################################################################################
#Gene regression parameter estimates

#Summary table: load previous save (TRUE); re-calculate (FALSE)
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
  library(ggrepel) #for volcano plot gene labels
  library(magrittr)
  library(MAST) #renv::install("RGLab/MAST")
  
  #all objects in the environment
  load(file = "./Data/zlm_1b_mid_sav.RData")
  #just the summary (estimates) of zlm_glmer
  #> est <- readRDS(file = "./Data/zlm_1b_est_sav.rds")
  
} else {
  est <- summary(zlm_glmer, doLRT = TRUE)  
}

#####

#Save?
if (FALSE) {
  #all objects in the environment
  #> save.image(file = "./Data/zlm_1b_mid_sav.RData")
  #just the summary (estimates) of zlm_glmer
  #> saveRDS(est, file = "./Data/zlm_1b_est_sav.rds")
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
est_dt <- as_tibble(est[["datatable"]]) %>% 
  dplyr::arrange(primerid, contrast, component)
est_dt$component %<>% as.factor()
est_dt$contrast %<>% as.factor()

#Discrete, Continuous, S (combined: has z-score but not coef or CI), logFC
#Each gene has:
#intercept (C, D, S)
#age.group67-86yo (C, D, logFC, S)
#age.group67-86yo:days_post42 (C, D, logFC, S)
#days_post42 (C, D, logFC, S)
summary(est_dt$component)
summary(est_dt$contrast)

#Check genes with NaN or NA
est_alt <- est_dt %>% dplyr::filter(component != "S")
est_alt <- est_alt[is.na(est_alt$coef),]
levels(as.factor(est_alt$primerid))  
#Should match number of failures here:
summary(zlm_glmer@converged)

################################################################################
#ENSG numerical order like in the estimates table
#For matching ENSG codes to the other gene codes
raw2_gene %<>% dplyr::arrange(primerid) 

#From the test result objects, extract the logFC part of the genes
#make p-value from z-score, make negative log 10 of p-value

#...day (young unvaccinated to young vaccinated)
# contrast == "days_post42"
vacc_tb_f <- as_tibble(est[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
vacc_tb_f$fdr <- p.adjust(vacc_tb_f$p_val, method = "BH")
vacc_tb_f %<>% relocate("fdr", .after = "neg_log_10_p")

##########

#...age (young unvaccinated to old unvaccinated)
# contrast == "age.group67-86yo"
age_tb_f <- as_tibble(est[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "age.group67-86yo") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
age_tb_f$fdr <- p.adjust(age_tb_f$p_val, method = "BH")
age_tb_f %<>% relocate("fdr", .after = "neg_log_10_p")

##########

#...interaction (young vaccinated to old vaccinated)
#log-fold change
inter_tb_f <- as_tibble(est[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "logFC" & contrast == "age.group67-86yo:days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
inter_tb_f$fdr <- p.adjust(inter_tb_f$p_val, method = "BH")
inter_tb_f %<>% relocate("fdr", .after = "neg_log_10_p")

#discrete
inter_tb_d <- as_tibble(est[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "D" & contrast == "age.group67-86yo:days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
inter_tb_d$fdr <- p.adjust(inter_tb_d$p_val, method = "BH")
inter_tb_d <- inter_tb_d %>% 
  relocate(c("ci.lo", "coef", "ci.hi", "p_val", "fdr"), .after = "symbol") %>% 
  arrange(p_val)
# inter_tb_d$ci.lo <- signif(inter_tb_d$ci.lo, 3)
# inter_tb_d$coef <- signif(inter_tb_d$coef, 3)
# inter_tb_d$ci.hi <- signif(inter_tb_d$ci.hi, 3)
# inter_tb_d$p_val <- signif(inter_tb_d$p_val, 2)
# inter_tb_d$fdr <- signif(inter_tb_d$fdr, 3)
inter_tb_d$EntrezGene <- as.character(inter_tb_d$EntrezGene)

#continuous
inter_tb_c <- as_tibble(est[["datatable"]]) %>% 
  arrange(primerid) %>% 
  dplyr::filter(component == "C" & contrast == "age.group67-86yo:days_post42") %>% 
  mutate(p_val = norm_2sp(z)) %>% 
  mutate(neg_log_10_p = -log10(p_val)) %>% 
  bind_cols(raw2_gene[,c("EntrezGene","symbol")]) %>% 
  dplyr::relocate(c("EntrezGene","symbol"), .after = "primerid") 
inter_tb_c$fdr <- p.adjust(inter_tb_c$p_val, method = "BH")
inter_tb_c <- inter_tb_c %>% 
  relocate(c("ci.lo", "coef", "ci.hi", "p_val", "fdr"), .after = "symbol") %>% 
  arrange(p_val)
# inter_tb_c$ci.lo <- signif(inter_tb_c$ci.lo, 3)
# inter_tb_c$coef <- signif(inter_tb_c$coef, 3)
# inter_tb_c$ci.hi <- signif(inter_tb_c$ci.hi, 3)
# inter_tb_c$p_val <- signif(inter_tb_c$p_val, 2)
# inter_tb_c$fdr <- signif(inter_tb_c$fdr, 3)
inter_tb_c$EntrezGene <- as.character(inter_tb_c$EntrezGene)

################################################################################
if (FALSE) {
  write_csv(inter_tb_f, file = "./Data/inter_tb_f.csv")
  write_csv(inter_tb_d, file = "./Data/inter_tb_d.csv")
  write_csv(inter_tb_c, file = "./Data/inter_tb_c.csv")
}

################################################################################
#Vaccine
#p-value hist
ggplot(data = vacc_tb_f, aes(x=p_val)) +
  geom_histogram(binwidth = 0.025, boundary = 0, na.rm = TRUE, 
    fill = "black", colour = "white", size = 0.2) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.20), minor_breaks = seq(0, 1, by = 0.025)) +
  theme_bw() +
  xlab("p-values for Vaccine log Fold-Changes") +
  ylab("Frequency")
  #ggtitle("vacc hist")

#Volcano Plot 
ggplot(data = vacc_tb_f, aes(x=coef, y=neg_log_10_p)) +
  geom_point(na.rm = TRUE, size = 1.5, alpha = 0.6, shape = 16) +
  geom_label_repel(data = vacc_tb_f %>% filter(fdr < 0.005),
    aes(label = symbol), 
    color = "red", alpha = 0.8, nudge_y = 1, force = 1, max.overlaps = 50
  ) + 
  coord_cartesian(xlim = c(-5,5), ylim = c(0, 11)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), minor_breaks = seq(-6, 6, by = 1)) +
  scale_y_continuous(breaks = seq(0, 11, by = 1), minor_breaks = NULL) +
  theme_bw() +
  xlab("log2(Fold-Change) for Vaccination") +
  ylab("-log10(p-value)") #+
  #ggtitle("vacc volc")


##########
#Age
#p-value hist
ggplot(data = age_tb_f, aes(x=p_val)) +
  geom_histogram(binwidth = 0.025, boundary = 0, na.rm = TRUE, 
    fill = "black", colour = "white", size = 0.2) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.20), minor_breaks = seq(0, 1, by = 0.025)) +
  theme_bw() +
  xlab("p-values for Age log Fold-Changes") +
  ylab("Frequency") #+
  #ggtitle("age hist")

#Volcano Plot 
ggplot(data = age_tb_f, aes(x=coef, y=neg_log_10_p)) +
  geom_point(na.rm = TRUE, size = 1.5, alpha = 0.6, shape = 16) +
  geom_label_repel(data = age_tb_f %>% filter(fdr < 0.05),
    aes(label = symbol), 
    color = "red", alpha = 0.8, nudge_y = 1, force = 1, max.overlaps = 30
  ) + 
  coord_cartesian(xlim = c(-5,5), ylim = c(0, 11)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), minor_breaks = seq(-6, 6, by = 1)) +
  scale_y_continuous(breaks = seq(0, 11, by = 1), minor_breaks = NULL) +
  theme_bw() +
  xlab("log2(Fold-Change) for Age") +
  ylab("-log10(p-value)") #+
  #ggtitle("age volc")

##########

#Interaction
#p-value hist
ggplot(data = inter_tb_f, aes(x=p_val)) +
  geom_histogram(binwidth = 0.025, boundary = 0, na.rm = TRUE, 
    fill = "black", colour = "white", size = 0.2) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.20), minor_breaks = seq(0, 1, by = 0.025)) +
  theme_bw() +
  xlab("p-values for Interaction log Fold-Changes") +
  ylab("Frequency") #+
  #ggtitle(intervolc)

#Volcano Plot log fold change
ggplot(data = inter_tb_f, aes(x=coef, y=neg_log_10_p)) +
  geom_point(na.rm = TRUE, size = 1.5, alpha = 0.6, shape = 16) +
  geom_label_repel(data = inter_tb_f %>% filter(fdr < 0.03),
    aes(label = symbol), 
    color = "red", alpha = 0.8, nudge_y = 2, nudge_x = 0, force = 2, force_pull = 3, max.overlaps = 30
  ) + 
  coord_cartesian(xlim = c(-5,5), ylim = c(0, 11)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), minor_breaks = seq(-6, 6, by = 1)) +
  scale_y_continuous(breaks = seq(0, 11, by = 1), minor_breaks = NULL) +
  theme_bw() +
  xlab("log2(Fold-Change) for interaction") +
  ylab("-log10(p-value)") #+
  #ggtitle("inter volc log fold change")

#####

#Volcano Plot discrete
ggplot(data = inter_tb_d, aes(x=coef, y=neg_log_10_p)) +
  geom_point(na.rm = TRUE, size = 1.5, alpha = 0.6, shape = 16) +
  geom_label_repel(data = inter_tb_d %>% filter(fdr < 0.05),
    aes(label = symbol), 
    color = "black", alpha = 0.8, nudge_y = 2, nudge_x = 0, force = 2, force_pull = 3, max.overlaps = 30
  ) + 
  coord_cartesian(xlim = c(-5,5), ylim = c(0, 11)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), minor_breaks = seq(-6, 6, by = 1)) +
  scale_y_continuous(breaks = seq(0, 11, by = 1), minor_breaks = NULL) +
  theme_bw() +
  xlab("log(Odds Ratio) for Interaction") +
  ylab("-log10(p-value)") #+
  #ggtitle("inter volc disc")

#Volcano Plot continuous
ggplot(data = inter_tb_c, aes(x=coef, y=neg_log_10_p)) +
  geom_point(na.rm = TRUE, size = 1.5, alpha = 0.6, shape = 16) +
  geom_label_repel(data = inter_tb_c %>% filter(fdr < 0.05),
    aes(label = symbol), 
    color = "red", alpha = 0.8, nudge_y = 2, nudge_x = 0, force = 2, force_pull = 3, max.overlaps = 30
  ) + 
  coord_cartesian(xlim = c(-5,5), ylim = c(0, 11)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), minor_breaks = seq(-6, 6, by = 1)) +
  scale_y_continuous(breaks = seq(0, 11, by = 1), minor_breaks = NULL) +
  theme_bw() +
  xlab("Mean Difference in et for Interaction") +
  ylab("-log10(p-value)") #+
  #ggtitle("inter cont")

####################

#discrete and continuous twoway
#make a temporary table of 

temp_d <- inter_tb_d %>% 
  arrange(primerid) %>% 
  select(c(primerid, EntrezGene, symbol, coef)) %>% 
  dplyr::rename(coef_d = coef)

temp_c <- inter_tb_c %>% 
  arrange(primerid) %>% 
  select(c(primerid, EntrezGene, symbol, coef)) %>% 
  dplyr::rename(coef_c = coef)

temp_f <- inter_tb_f %>% 
  arrange(primerid) %>% 
  select(c(primerid, EntrezGene, symbol, coef, p_val, fdr)) %>% 
  dplyr::rename(coef_f = coef, p_val_f = p_val, fdr_f = fdr)

temp_d_vs_c <- temp_d %>% 
  bind_cols("coef_c" = temp_c$coef_c) %>% 
  bind_cols("coef_f" = temp_f$coef_f) %>% 
  bind_cols("p_val_f" = temp_f$p_val_f) %>% 
  bind_cols("fdr_f" = temp_f$fdr_f)
  
ggplot(data = temp_d_vs_c) +
  geom_point(aes(x=coef_d, y=coef_c), na.rm = TRUE, size = 1.5, alpha = 0.6, shape = 16) +
  xlab("log odds ratio \n (old age, vaccine (after): interaction)") +
  ylab("mean difference in et \n (old age, vaccine (after): interaction)") +
  geom_label_repel(data = temp_d_vs_c %>% filter(fdr_f < 0.03),
    aes(x=coef_d, y=coef_c, label = symbol), 
    color = "red", alpha = 0.8, nudge_y = 0, nudge_x = -2, force = 2, force_pull = 3, max.overlaps = 30
  ) +
  theme_bw() #+
  #ggtitle("inter volc")

################################################################################
#Fitted baseline log odds of passing hurdle (young, day0)
b1 <- tibble("coefD_base" = as.vector(zlm_glmer@coefD[,1]))
ggplot(data = b1, aes(x=coefD_base)) +
  geom_histogram(boundary = 0, binwidth = 0.25, na.rm = TRUE) +
  scale_x_continuous(breaks = seq(-10, 5, by = 5), minor_breaks = seq(-10, 5, by = 1)) +
  theme_bw()

################################################################################
#Test gene groups for the interaction term
#Sort in p-value order
inter_tb_f_psort <- inter_tb_f %>% dplyr::arrange(p_val) 
inter_tb_d_psort <- inter_tb_d %>% dplyr::arrange(p_val) 
inter_tb_c_psort <- inter_tb_c %>% dplyr::arrange(p_val) 

top_f <- inter_tb_f_psort %>% dplyr::filter(fdr < 0.05)
top_d <- inter_tb_d_psort %>% dplyr::filter(fdr < 0.05) #empty
top_c <- inter_tb_c_psort %>% dplyr::filter(fdr < 0.05)

if (FALSE) {
  write_csv(top_f, file = "./Data/2btop_f.csv")
  write_csv(inter_tb_d_psort, file = "./Data/2btop_dz.csv")
  write_csv(top_c, file = "./Data/2btop_c.csv")
}

#renv::install("bioc::GO.db")
library(limma)

top_f_go <- tibble(goana(top_f$EntrezGene, species = "Hs")) %>% dplyr::arrange(P.DE)
top_f_go_cc <- top_f_go %>% dplyr::filter(Ont == "CC")
top_f_go_bp <- top_f_go %>% dplyr::filter(Ont == "BP")
top_f_go_mf <- top_f_go %>% dplyr::filter(Ont == "MF")

top_f_kg <- tibble(kegga(top_f$EntrezGene, species = "Hs")) %>% dplyr::arrange(P.DE)

if (FALSE) {
  write_csv(top_f_go_cc, file = "./Data/2btop_f_go_cc.csv")
  write_csv(top_f_go_bp, file = "./Data/2btop_f_go_bp.csv")
  write_csv(top_f_go_mf, file = "./Data/2btop_f_go_mf.csv")
  write_csv(top_f_kg, file = "./Data/2btop_f_kg.csv")
}
