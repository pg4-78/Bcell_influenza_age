################################################################################
#Script-file: diagnostic2.R
#Description: histogram and cdf to show excess zeros
################################################################################

#immediately after setup1a, setup1b
#or alternatively, load save
if(TRUE) {
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
  load(file = "./Data/setup1_b3_fltr_sav.RData")
}

################################################################################
#Check number of samples per individual
table(raw2_samp$indiv, raw2_samp$days_post)

#Check which individuals belong to each age.group
table(raw2_samp$indiv, raw2_samp$age.group)

################################################################################
#Check et
#The gene filter cutoff points chosen in setup1b are 
#... 0.10 for minimum proportion of samples expressing the gene 
#... 0.6 for minimum average et of samples for the gene
dim(raw2_mx_et)
dim(sca_l)

#Put row numbers and inspect
raw2_gene <- bind_cols(raw2_gene, "num" = 1:dim(raw2_gene)[1])

#Try specific genes with moderate amounts of et
#gene in row 5435 has avg et = 1.0001727
#gene in row 5077 has avg et = 2.000014
#gene in row 2015 has avg et = 3.000463
#gene in row 1405 has avg et = 4.0096832

#genes of interest
#choose a specific gene (TRUE), use all genes (FALSE)


if (TRUE) {
  g <- 1405
} else {
  g <- 1:dim(raw2_gene)[1]
}

#samples of interest
#choose a specific individual and day (TRUE), use all samples (FALSE)
if (FALSE) {
  s <- rep(NA, dim(raw2_samp)[1])
  for (i in 1:dim(raw2_samp)[1]) {
    if (raw2_samp$indiv[i]=="543P" & raw2_samp$days_post[i]=="42") {
      s[i] <- i
    }
  }
  s <- s[!is.na(s)]
} else {
  s <- 1:dim(raw2_samp)[1]
}

#check et histogram
ggplot(data = tibble("et" = as.vector(raw2_mx_et[g,s])), aes(x=et)) +
  geom_histogram(boundary = 0.000001, binwidth = 0.5, fill="black", colour="white", size=0.2) +
  scale_x_continuous(breaks = seq(0, 20, by = 5), minor_breaks = seq(0, 20, by = 0.5)) +
  geom_vline(xintercept=0, colour = "red", size = 1) +
  geom_hline(yintercept=0, colour = "black", size = 0.5) +
  coord_cartesian(x = c(0, 12), y = c(0, 650)) + 
  theme_bw() +
  xlab("et (Transformed Count)") +
  ylab("Frequency")

#check et cdf
ggplot(data = tibble("et" = round(as.vector(raw2_mx_et[g,s]), digits = 3)), aes(x=et)) +
  stat_ecdf() +
  scale_x_continuous(breaks = seq(0, 20, by=5), minor_breaks = seq(0, 20, by = 0.5)) +
  scale_y_continuous(minor_breaks = seq(0, 1, by=0.1), breaks = seq(0, 1, by=0.2)) +
  coord_cartesian(x = c(0, 12), y = c(0, 1)) + 
  theme_bw() +
  xlab("et (Transformed Count)") +
  ylab("Cumulative Distribution Function")

####################
#Check qq plot
a1 <- tibble("et" = as.vector(raw2_mx_et[g,s]))
a2 <- a1 %>% filter (et > 0)

ggplot(data = a1, aes(sample=et)) +
  geom_qq() +
  geom_qq_line(col = "red") +
  theme_bw()
  
ggplot(data = a2, aes(sample=et)) +
  geom_qq() +
  geom_qq_line(col = "red") + 
  theme_bw()
