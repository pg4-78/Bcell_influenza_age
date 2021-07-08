################################################################################
#Script-file: setup1b.R
#Description: data collation PART B
################################################################################

#Directly after running setup1a and answering prompt (database selection)
#or alternatively
if (TRUE) {
  #Clear
  #...variables
  rm(list=ls())
  #...console
  cat("\014\n")
  #...graphs
  tryCatch(dev.off(), error = function(e) {NULL})
  dev.new()
  
  load(file = "./Data/setup1_a3_sav.RData") #> 
  
  library(tidyverse)
  library(magrittr)
  #library(org.Hs.eg.db) 
  library(EDASeq)
  library(GEOquery)
  library(SingleCellExperiment)
  library(MAST) #renv::install("RGLab/MAST@bb1e928c2fc7bcd2fdbf19472e43bf6848e9c4dc")
  #Keep "RSQLite@2.2.5"
}

################################################################################
#Attach length and GC content to gene information table
gene_info %<>%  as_tibble()
raw_gene %<>% bind_cols(gene_info)

sum(!is.na(raw_gene$length)) #Strands with lengths
sum(is.na(raw_gene$length)) #Strands without lengths


#####
#Column totals before gene filter
raw_mx_bef_fltr <- raw %>% dplyr::select(-"Ensembl")

tot_v_bef_fltr <- rep(0, dim(raw_mx_bef_fltr)[2])

for (col in 1:dim(raw_mx_bef_fltr)[2]) {
  tot_v_bef_fltr[col] <- sum(raw_mx_bef_fltr[,col])
}

########################################
#raw2 eventually gets split into count matrix (/transformed) and gene details
#####
#Filter to include only those with entries for length
#Normalise counts by gene length and total counts in each sample

raw2 <-  bind_cols(raw_gene, (raw %>% dplyr::select(-"Ensembl"))) %>% 
  filter(!is.na(length)) %>% 
  dplyr::mutate("klength" = length/1000) %>% 
  relocate(klength, .after = length)

#Part with just the counts, calculate et from counts
raw2_mx_ct <- as.matrix(raw2 %>% dplyr::select(-c("Ensembl", "EntrezGene", "symbol", "length", "klength", "gc")))
raw2_mx_et <- raw2_mx_ct
raw2_klength <- as.vector(raw2$klength)

#####
#Column totals after gene filter
tot_v_aft_fltr <- rep(0, dim(raw2_mx_ct)[2])

for (col in 1:dim(raw2_mx_ct)[2]) {
  tot_v_aft_fltr[col] <- sum(raw2_mx_ct[,col])
}

#Relationship of cell (column) count total after against before
#filter for only genes with no length
sum(tot_v_aft_fltr <= tot_v_bef_fltr)

ggplot(data = tibble("bef" = tot_v_bef_fltr, "aft" = tot_v_aft_fltr)) +
  geom_point(aes(x=bef,y=aft), size = 1) +
  theme_bw() +
  xlab("Cell (Column) Count Total Before") +
  ylab("Cell (Column) Count Total After") +
  geom_abline(slope=1, intercept=0, colour="red", size=0.7)

#####
#Normalise 1: Divide each row by gene length (KILObases)
for (row in 1:dim(raw2_mx_et)[1]) {
  if (sum(raw2_mx_et[row,]) > 0) {
    raw2_mx_et[row,] <- raw2_mx_et[row,] / raw2_klength[row]
  }
}

#Normalise 2: Each column (sample) gets divided by (column count total in MILLIONS)
for (col in 1:dim(raw2_mx_et)[2]) {
  raw2_mx_et[,col] <- raw2_mx_et[,col] / (sum(raw2_mx_et[,col]) / 10^6)
}

#The matrix is now in the (T)ranscripts (P)er (M)illion scale
sum(raw2_mx_et[,1])
sum(raw2_mx_et[,2])

#MAST uses et = log2(TPM + 1)
raw2_mx_et <- log2(raw2_mx_et+1)

########################################
#Feature (gene) data with as many rows as expression matrix rows
#...only genes with length available
raw2_gene <- raw2 %>% dplyr::select(c("Ensembl", "EntrezGene", "symbol", "length", "klength", "gc")) %>% 
  dplyr::rename(primerid = Ensembl)

########################################
#For ease of filtering, get:
#fraction of samples with the gene expressed
#average expression level across samples (including the zeros)

raw2_gene %<>% 
  mutate(frac_active = 0) %>% 
  mutate(avg_et = 0)

#Check dimension match b/w gene info and counts
dim(raw2_gene)[1] == dim(raw2_mx_et)[1]
dim(raw2_gene)[1] == dim(raw2_mx_ct)[1]

for (i in 1:(dim(raw2_gene)[1])) { 
  raw2_gene$frac_active[i] <- sum(raw2_mx_et[i, ]>0.000001)/dim(raw2_mx_et)[2]
  raw2_gene$avg_et[i] <- sum(raw2_mx_et[i, ])/dim(raw2_mx_et)[2]
}

quantile(raw2_mx_et, 0.91)
summary(raw2_gene$frac_active)
summary(raw2_gene$avg_et)

#Step to filter by these levels towards the end

########################################
#Cell/ sample data (as many rows as expression matrix columns)
n_gene <- dim(raw2_mx_et)[1]
n_samp <- dim(raw2_mx_et)[2]

#Get all the experiment data
gse <- getGEO("GSE167823", GSEMatrix = FALSE)
length(gse)
names(gse)

#Extract individual samples 
names(GSMList(gse))
#Try the first individual
GSMList(gse)[[1]]

#Extract the individual conditions/ characteristics: 
#...(eg. first sample) before looping over all samples
gse@gsms[["GSM5112843"]]@header[["characteristics_ch1"]]
str_split(gse@gsms[["GSM5112843"]]@header[["title"]], pattern = "_")

#Loop through and get all of these into a matrix
#...1 treatment; 2 time; 3 age; 4 age.group;
#...5 cell type; 6 cell type_sorted

raw2_samp <- tibble(
  "treatment" = rep("0", n_samp),
  "days_post" = rep("0", n_samp), #Time, days_post
  "age" = rep("0", n_samp),
  "age.group" = rep("0", n_samp),
  "cell type" = rep("0", n_samp),
  "cell type_sorted" = rep("0", n_samp),
  #
  "indiv" = rep("0", n_samp),
  "cDNA" = rep("0", n_samp),
  "well" = rep("0", n_samp),
  "wellKey" = colnames(raw)[-1],
  "CDR" = rep(0, n_samp)
)

for (i in 1:n_samp) {
  raw2_samp[i,"treatment"] <- word(gse@gsms[[i]]@header[["characteristics_ch1"]][[1]], start = 2, end = -1)
  raw2_samp[i,"days_post"] <- word(gse@gsms[[i]]@header[["characteristics_ch1"]][[2]], start = 3, end = -1) #Time, days_post
  raw2_samp[i,"age"] <- word(gse@gsms[[i]]@header[["characteristics_ch1"]][[3]], start = 2, end = -1)
  raw2_samp[i,"age.group"] <- word(gse@gsms[[i]]@header[["characteristics_ch1"]][[4]], start = 2, end = -1)
  raw2_samp[i,"cell type"] <- word(gse@gsms[[i]]@header[["characteristics_ch1"]][[5]], start = 3, end = -1)
  raw2_samp[i,"cell type_sorted"] <- word(gse@gsms[[i]]@header[["characteristics_ch1"]][[6]], start = 3, end = -1)

  temp <- str_split(gse@gsms[[i]]@header[["title"]], pattern = "_")
  raw2_samp[[i,"indiv"]] <- temp[[1]][1]
  raw2_samp[[i,"cDNA"]] <- temp[[1]][4]
  raw2_samp[[i,"well"]] <- temp[[1]][5]
}

#treatment 
raw2_samp$treatment %<>% forcats::as_factor()
raw2_samp$treatment %<>% plyr::revalue(c(`trivalent influenza vaccination` = "vaccine"))
summary(raw2_samp$treatment)

#days_post
raw2_samp$days_post %<>% forcats::as_factor() #Time, days_post
summary(raw2_samp$days_post) #Time, days_post

#age 
raw2_samp$age %<>% as.numeric()
summary(raw2_samp$age)
summary(forcats::as_factor(raw2_samp$age))

#age.group 
raw2_samp$age.group %<>% forcats::as_factor()
summary(raw2_samp$age.group)

#`cell type` ... constant across all samples
summary(forcats::as_factor(raw2_samp$`cell type`))

#`cell type_sorted` ... constant across all samples
summary(forcats::as_factor(raw2_samp$`cell type_sorted`))

#indiv
raw2_samp$indiv <- forcats::as_factor(raw2_samp$indiv)
summary(raw2_samp$indiv)

#cDNA
if (FALSE) {raw2_samp$cDNA %<>% as.character()} #RESET
raw2_samp$cDNA %<>% forcats::as_factor()
summary(raw2_samp$cDNA)

#well 
if (FALSE) {raw2_samp$well %<>% as.character()} #RESET
raw2_samp$well %<>% factor(levels = 
  c(
    "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",
    "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12",
    "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
    "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
    "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12",
    "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12",
    "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12",
    "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12"
  )
)
summary(raw2_samp$well)
levels(raw2_samp$well)

################################################################################
#CDR: proportion of genes above 0
#Loop through each sample (matrix column)
#Add it to the CDR column (as a sample/cell/column detail)
for (j in 1:dim(raw2_mx_et)[2]) {
  raw2_samp$CDR[j] = sum(raw2_mx_et[,j]!=0)/(dim(raw2_mx_et)[1])
}
summary(raw2_samp$CDR)

################################################################################
#For testing: alternate coding
#make an age factor variable recoded as 1 young, 0 old
#make a vaccine factor variable recoded as 1 before, 0 after
raw2_samp %<>% mutate("before" = ifelse(days_post==0,1,0))
raw2_samp$before %<>% forcats::as_factor()

raw2_samp %<>% mutate("young" = ifelse(age.group=="22-36yo",1,0))
raw2_samp$young %<>% forcats::as_factor()

#Double check the alternative coding
table(raw2_samp$days_post, raw2_samp$before)
table(raw2_samp$age.group, raw2_samp$young)

################################################################################
#Filter genes: cut low fraction active and low average et

#Before
dim(raw2_mx_et)[1]
dim(raw2_gene)[1] == dim(raw2_mx_et)[1]
dim(raw2_gene)[1] == dim(raw2_mx_ct)[1]
#
mean(raw2_gene$frac_active)
sd(raw2_gene$frac_active)
quantile(raw2_gene$frac_active, c(0, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 1))
#
mean(raw2_gene$avg_et)
sd(raw2_gene$avg_et)
quantile(raw2_gene$avg_et, c(0, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 1))

#Cutoff points
min_frac_active <- 0.10 #quantile(raw2_gene$frac_active, 0)
max_frac_active <- quantile(raw2_gene$frac_active, 1)
min_avg_et <-  0.6 #quantile(raw2_gene$avg_et, 0)
max_avg_et <- quantile(raw2_gene$avg_et, 1)

include_gene <- (
  (raw2_gene$frac_active >= min_frac_active) & 
  (raw2_gene$frac_active <= max_frac_active) & 
  (raw2_gene$avg_et >= min_avg_et) &
  (raw2_gene$avg_et <= max_avg_et)
)

#(Transformed) count matrix
raw2_mx_et <- raw2_mx_et[include_gene==TRUE,]

#Original count matrix
raw2_mx_ct <- raw2_mx_ct[include_gene==TRUE,]

#Gene list 
raw2_gene <- raw2_gene[include_gene==TRUE,]

#After
dim(raw2_mx_et)[1]
dim(raw2_gene)[1] == dim(raw2_mx_et)[1]
dim(raw2_gene)[1] == dim(raw2_mx_ct)[1]
#
mean(raw2_gene$frac_active)
sd(raw2_gene$frac_active)
quantile(raw2_gene$frac_active, c(0, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 1))
#
mean(raw2_gene$avg_et)
sd(raw2_gene$avg_et)
quantile(raw2_gene$avg_et, c(0, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 1))

################################################################################
#Checks 

#If the column variable uniquely identifies individuals,
#min/max ages should be the same (or differ by 1 year)
for (i in levels(raw2_samp$indiv)) {
  print(i)
  print(summary(raw2_samp[raw2_samp$indiv == i, "age"]))
  print("----------")
}

#The same cDNA can be found in multiple individuals
for (i in levels(raw2_samp$cDNA)) {
  print(i)
  print(summary(raw2_samp[raw2_samp$cDNA == i, "age"]))
  print("----------")
}

#Check for duplicates
if (FALSE) {
  View(raw2_gene[duplicated(raw2_gene$EntrezGene)==1,])
}

#For each individual, check that there are cells present for day 0 and 42
#...All the combinations of individual and day have samples
for (i in levels(raw2_samp$indiv)){
  print(i)
  print(summary(forcats::as_factor(raw2_samp[raw2_samp$indiv==i, "days_post"]))) #Time, days_post
  print("----------")
}

#Age group by individual
table(raw2_samp$age.group, raw2_samp$indiv)

#Age group by vaccine status
table(raw2_samp$age.group, raw2_samp$days_post)

#Age by individual
for (i in levels(raw2_samp$indiv)){
  print(i)
  print(summary(forcats::as_factor(raw2_samp[raw2_samp$indiv==i, "age"])))
  print("----------")
}

################################################################################
#Gather all the required parts to form the SingleCellAssay object
#[primerid, wellKey] are just for subsetting

#Longitudinal version
sca_l <- FromMatrix(
  exprsArray = raw2_mx_et, 
  cData = as.data.frame(raw2_samp), 
  fData = as.data.frame(raw2_gene)
)

#assays(sca_l, withDimnames=FALSE)$counts <- raw2_mx_ct
assayNames(sca_l)

if(FALSE) {
  View(sca_l@assays@data@listData[["et"]], title = "sca_l_et")
  View(sca_l@assays@data@listData[["counts"]], title = "sca_l_ct")
}
#Cross-sectional version include only observations at day 42
#...exclude those from day 0
sca_d42 <- BiocGenerics::subset(sca_l, days_post==42) #Time, days_post
sca_d00 <- BiocGenerics::subset(sca_l, days_post==0) #Time, days_post
#... double check that they add to original sample size (952)
dim(sca_d42)[2] + dim(sca_d00)[2]

################################################################################
#Clean up
rm(n_gene, n_samp, row, col, i, j, temp, m, raw2_klength)

#Save?
if (FALSE) {
  save.image(file = "./Data/setup1_b3_fltr_sav.RData")
}

if (FALSE) {
  #updated with reverse coding included
  #(young 1, old 0); (before 1, after 0)
  sca_l_rev <- sca_l
  rm(sca_l)
  save.image(file = "./Data/setup1_b3_fltr_rev_sav.RData")
}

#setup1_b1_sav: before MAST update
#setup1_b2_sav: with MAST update
#setup1_b3_sav: also with BSgenome, TxDb update

#no low expr filter yet: goes to zlm_test1
#...to figure out what conditions make the model fail to converge
#file = "./Data/setup1_b3_sav.RData" 

#low expr filter active: goes to zlm1a
#using the options: 
#min_frac_active: 0.10 #min_avg_et: 0.6 
#file = "./Data/setup1_b3_fltr_sav.RData"

#also with reverse coding for testing
#(young 1, old 0); (before 1, after 0)
#file = "./Data/setup1_b3_fltr_rev_sav.RData"