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
  
  #Genes filtered by frac_active and avg_et
  load(file = "./Data/setup1_b3_fltr_sav.RData")
}

#Use all cores for faster computation?
if (TRUE) {
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
# rm(list=setdiff(ls(), c("list", "to", "keep")))
rm(gse, raw, raw_gene, egENSEMBL_tb, egSymbol_tb, gene_info)
rm(core_lim, opt_temp)

#Try timing first 100 genes
if (TRUE) {
  time_s <- Sys.time()
  zlm_glmer <- zlm(
    ~age.group*days_post + (1|indiv),
    sca = sca_l[1:100,],
    method = "glmer", ebayes = FALSE, parallel = TRUE
  )
  time_f <- Sys.time()
  print(time_f - time_s)
  rm(time_s, time_f)
}

################################################################################
#Checks
dim(raw2_mx_et)
dim(sca_l)

if (FALSE) {
  summary(sca_d00@colData@listData[["age.group"]])
  summary(sca_d42@colData@listData[["age.group"]])
  summary(sca_l@colData@listData[["age.group"]])
  #####
  summary(sca_d00@colData@listData[["indiv"]])
  summary(sca_d42@colData@listData[["indiv"]])
  summary(sca_l@colData@listData[["indiv"]])
}

################################################################################
#Using zlm SingleCellAssay glmer
#parallel to use multiple cores
zlm_glmer <- zlm(
  ~age.group*days_post + (1|indiv),
  sca = sca_l,
  method = "glmer", ebayes = FALSE, parallel = TRUE
)

#The previous function takes long 
#...Save?
if (FALSE) {
  save.image(file = "./Data/zlm1a_sav.RData")
}
