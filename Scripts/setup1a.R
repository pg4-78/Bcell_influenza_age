#raw: original dataset
#raw2: only genes with length data available
#SingleCellAssay objects (only genes with length data available):
#sca_l: longitudinal: both days 0 and 42
#sca_d42: cross-sectional: day 42 only
#sca_d00: cross-sectional: day 0 only
################################################################################

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
#library(org.Hs.eg.db) 
library(EDASeq)
library(GEOquery)
library(SingleCellExperiment)
library(MAST) #renv::install("RGLab/MAST@bb1e928c2fc7bcd2fdbf19472e43bf6848e9c4dc")
#Keep "RSQLite@2.2.5"

################################################################################
#Data sub-folder inside the project
raw <- readr::read_csv("./Data/GEO_supporting_processed_data_file_raw_count_matrix.csv")
colnames(raw)[1] = "Ensembl"

################################################################################
#The MAST compatible SingleCellAssay object 

########################################
#Feature (gene) data with as many rows as expression matrix rows
#Use databases to get Entrez and symbols corresponding to Ensembl 
#Also find the corresponding gene lengths
raw_gene <- tibble("Ensembl" = raw$Ensembl)

egENSEMBL_tb <- toTable(org.Hs.egENSEMBL)
m <- match(raw_gene$Ensembl, egENSEMBL_tb$ensembl_id)
raw_gene$EntrezGene <- egENSEMBL_tb$gene_id[m]

egSymbol_tb <- toTable(org.Hs.egALIAS2EG)
m <- match(raw_gene$EntrezGene, egSymbol_tb$gene_id)
raw_gene$symbol <- egSymbol_tb$alias_symbol[m]

#Require the gene length
#... options: mode = c("biomart", "org.db")
#... choose: "TxDb.Hsapiens.UCSC.hg19.knownGene", "BSgenome.Hsapiens.UCSC.hg19"
#Keep "RSQLite@2.2.5"
#> save.image(file = "./Data/setup1_a3_sav.RData") 
#...after finishing
#Continue to setup1b
gene_info <- getGeneLengthAndGCContent(raw_gene$Ensembl, org = "hsa", mode = "org.db")
