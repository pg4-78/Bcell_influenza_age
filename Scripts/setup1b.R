#Directly after running setup1a and answering prompt (database selection)

gene_info %<>%  as_tibble()
raw_gene %<>% bind_cols(gene_info)

sum(!is.na(raw_gene$length)) #Strands with lengths
sum(is.na(raw_gene$length)) #Strands without lengths

########################################
#Filter to include only those with entries for length
#Normalise counts by gene length and total counts in each sample
raw2 <-  bind_cols(raw_gene, (raw %>% dplyr::select(-"Ensembl"))) %>% 
  filter(!is.na(length)) %>% 
  dplyr::mutate("klength" = length/1000) %>% 
  relocate(klength, .after = length)
  
raw2_mx_ct <- as.matrix(raw2 %>% dplyr::select(-c("Ensembl", "EntrezGene", "symbol", "length", "klength", "gc")))
raw2_mx_et <- raw2_mx_ct
raw2_klength <- as.vector(raw2$klength)

#Normalise 1: Divide each row by gene length (KILObases) #
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
#Each column gets divided by the total 

########################################
#Feature (gene) data with as many rows as expression matrix rows
#...only genes with length available
raw2_gene <- raw2 %>% dplyr::select(c("Ensembl", "EntrezGene", "symbol", "length", "klength", "gc")) %>% 
  dplyr::rename(primerid = Ensembl)

########################################
#Cell data (as many rows as expression matrix columns)
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
  "days_post" = rep("0", n_samp),
  "age" = rep("0", n_samp),
  "age.group" = rep("0", n_samp),
  "cell type" = rep("0", n_samp),
  "cell type_sorted" = rep("0", n_samp),
  #
  "indiv" = rep("0", n_samp),
  "cDNA" = rep("0", n_samp),
  "well" = rep("0", n_samp),
  "wellKey" = colnames(raw)[-1]
)

for (i in 1:n_samp) {
  raw2_samp[i,"treatment"] <- word(gse@gsms[[i]]@header[["characteristics_ch1"]][[1]], start = 2, end = -1)
  raw2_samp[i,"days_post"] <- word(gse@gsms[[i]]@header[["characteristics_ch1"]][[2]], start = 3, end = -1)
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
raw2_samp$days_post %<>% forcats::as_factor()
summary(raw2_samp$days_post)

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
raw2_samp$indiv %<>% forcats::as_factor()
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

####################
#______Checks

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
#...All the combinations of individual:day have samples
for (i in levels(raw2_samp$indiv)){
  print(i)
  print(summary(forcats::as_factor(raw2_samp[raw2_samp$indiv==i, "days_post"])))
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

assays(sca_l, withDimnames=FALSE)$counts <- raw2_mx_ct
assayNames(sca_l)

if(FALSE) {
  View(sca_l@assays@data@listData[["et"]], title = "sca_l_et")
  View(sca_l@assays@data@listData[["counts"]], title = "sca_l_ct")
}
#Cross-sectional version include only observations at day 42
#...exclude those from day 0
sca_d42 <- BiocGenerics::subset(sca_l, days_post==42)
sca_d00 <- BiocGenerics::subset(sca_l, days_post==0)
#... double check that they add to original sample size (952)
dim(sca_d42)[2] + dim(sca_d00)[2]

################################################################################
#Clean up
rm(n_gene, n_samp, row, col, i, temp, m, raw2_klength)

#Save?
if (FALSE) {
  save.image(file = "./Data/setup1_a_sav.RData")
}
