################################################################################
#Script-file: diagnostic1.R
#Description: PCA plots
################################################################################

#Immediately after setup1a and setup1b 
#or alternatively, load a previous save
if (TRUE) {
  #Clear
  #...variables
  rm(list=ls())
  #...console
  cat("\014\n")
  #...graphs
  tryCatch(dev.off(), error = function(e) {NULL})
  dev.new()
  
  load(file = "./Data/setup1_b3_fltr_sav.RData")
  #file = "./Data/setup1_b3_sav.RData"
  #file = "./Data/setup1_b3_fltr_sav.RData"
}

library(magrittr)
library(tidyverse)
library(EDASeq)

################################################################################
#using the calculation procedure of findMethods(plotPCA)
#but display the coordinates with ggplot instead

coord_PCA_fn <- function (object, k = 2, labels = TRUE, isLog = TRUE, ...) {
  if (!isLog) {
      Y <- apply(log(object + 1), 1, function(y) scale(y, 
          center = TRUE, scale = FALSE))
  }
  else {
      Y <- apply(object, 1, function(y) scale(y, center = TRUE, scale = FALSE))
  }
  s <- svd(Y)
  percent <- s$d^2/sum(s$d^2) * 100
  labs <- sapply(
    seq_along(percent), 
    function(i) {
      paste("PC ", i, " (", round(percent[i], 2), "%)", sep = "")
    }
  )
  return(c("s" = s, "labs" = labs))
}
################################################################################
#Select samples
#Try, for example, extracting all samples from one individual: eg.
#> pca_samp_select <- raw2_samp$indiv=="501T"

#age.group = c("22-36yo", "67-86yo")
#days_post = c("0", "42")
#indiv = c(
#"501T", "520P"*o2, "526W"*o3, "536G", "541M", 
#"543P"*o1, "544Q", "545R", "559G", "562K"*y3, 
#"568R", "594V", "602D", "622A", "627F"*y2, 
#"637R", "643Y", "652H", "660R", "665X"*y1
#)

#         501T 520P 526W 536G 541M
# 22-36yo   48    0    0    0    0
# 67-86yo    0   48   48   42   43
# 
#         543P 544Q 545R 559G 562K
# 22-36yo    0    0    0   54   59
# 67-86yo   72   48   45    0    0
# 
#         568R 594V 602D 622A 627F
# 22-36yo   48   59   34   48   64
# 67-86yo    0    0    0    0    0
# 
#         637R 643Y 652H 660R 665X
# 22-36yo   38    0    0    0   88
# 67-86yo    0   19    6   41    0

#alternatively, to continue without removing any samples
#pca_samp_select <- rep(TRUE, dim(raw2_samp)[1])

pca_samp_select <- raw2_samp$indiv=="526W"


##########
#new table with only selected samples
pca_samp <- raw2_samp[pca_samp_select==TRUE, ]
#matrix used in the PCA calculation
pca_mx_et <- raw2_mx_et[,pca_samp_select==TRUE]

dim(pca_samp)
dim(pca_mx_et)

#Calculate principal components
#"x" = s$u[, 1], "y" = s$u[, 2]
coord_PCA_calc <- coord_PCA_fn(pca_mx_et) #Uses isLog TRUE

#Make a table to help with ggplot
coord_PCA_tb <- tibble(
  "pc1" = coord_PCA_calc[["s.u"]][,1], 
  "pc2" = coord_PCA_calc[["s.u"]][,2],
  "pc3" = coord_PCA_calc[["s.u"]][,3],
  "age.group" = pca_samp$age.group, 
  "days_post" = pca_samp$days_post, 
  "label" = "0",
  "indiv" = pca_samp$indiv,
  "colour_grp" = "0"
)

################################################################################
#Choose four colours for colouring by age and day
#Label

(coord_PCA_tb$colour_grp[coord_PCA_tb$age.group=="22-36yo" & coord_PCA_tb$days_post=="0"] 
  <- "#F0A400")
(coord_PCA_tb$label[coord_PCA_tb$age.group=="22-36yo" & coord_PCA_tb$days_post=="0"] 
  <- "young, day 0")

#####
(coord_PCA_tb$colour_grp[coord_PCA_tb$age.group=="67-86yo" & coord_PCA_tb$days_post=="0"] 
  <- "#42BAFF")
(coord_PCA_tb$label[coord_PCA_tb$age.group=="67-86yo" & coord_PCA_tb$days_post=="0"] 
  <- "old, day 0")

#####
(coord_PCA_tb$colour_grp[coord_PCA_tb$age.group=="22-36yo" & coord_PCA_tb$days_post=="42"] 
  <- "#F00000")
(coord_PCA_tb$label[coord_PCA_tb$age.group=="22-36yo" & coord_PCA_tb$days_post=="42"] 
  <- "young, day 42")

#####
(coord_PCA_tb$colour_grp[coord_PCA_tb$age.group=="67-86yo" & coord_PCA_tb$days_post=="42"] 
  <- "#5A00A3")
(coord_PCA_tb$label[coord_PCA_tb$age.group=="67-86yo" & coord_PCA_tb$days_post=="42"] 
  <- "old, day 42")

####################

coord_PCA_tb %<>% dplyr::arrange(desc(days_post))

ggplot(data = coord_PCA_tb, aes(x=pc1, y=pc2)) +
  geom_point(aes(color = colour_grp), size=1) +
  scale_color_identity(
    "Group", 
    labels = coord_PCA_tb$label,
    breaks = coord_PCA_tb$colour_grp,
    guide = "legend"
  ) +
  xlab(coord_PCA_calc[["labs1"]]) +
  ylab(coord_PCA_calc[["labs2"]]) +
  theme_bw() #+
  #ggtitle("coloured by age and day")

#can try other dimension pairs eg (pc1, pc3); (pc2, pc3)
#by changing aes(x= pc_, y= pc_)
#and "labs_"

################################################################################
#Just two colours by day (regardless of age)
#Label

(coord_PCA_tb$colour_grp[coord_PCA_tb$age.group=="22-36yo" & coord_PCA_tb$days_post=="0"] 
  <- "#FF1000")
(coord_PCA_tb$label[coord_PCA_tb$age.group=="22-36yo" & coord_PCA_tb$days_post=="0"] 
  <- "young, day 0")

#####
(coord_PCA_tb$colour_grp[coord_PCA_tb$age.group=="67-86yo" & coord_PCA_tb$days_post=="0"] 
  <- "#FF1000")
(coord_PCA_tb$label[coord_PCA_tb$age.group=="67-86yo" & coord_PCA_tb$days_post=="0"] 
  <- "old, day 0")

#####
(coord_PCA_tb$colour_grp[coord_PCA_tb$age.group=="22-36yo" & coord_PCA_tb$days_post=="42"] 
  <- "#000000")
(coord_PCA_tb$label[coord_PCA_tb$age.group=="22-36yo" & coord_PCA_tb$days_post=="42"] 
  <- "young, day 42")

#####
(coord_PCA_tb$colour_grp[coord_PCA_tb$age.group=="67-86yo" & coord_PCA_tb$days_post=="42"] 
  <- "#000000")
(coord_PCA_tb$label[coord_PCA_tb$age.group=="67-86yo" & coord_PCA_tb$days_post=="42"] 
  <- "old, day 42")

####################

ggplot(data = coord_PCA_tb, aes(x=pc1, y=pc2)) +
  geom_point(aes(color = colour_grp), size=1) +
  scale_color_identity(
    "Group", 
    labels = coord_PCA_tb$label,
    breaks = coord_PCA_tb$colour_grp,
    guide = FALSE
  ) +
  xlab(coord_PCA_calc[["labs1"]]) +
  ylab(coord_PCA_calc[["labs2"]]) +
  theme_bw() #+
  #ggtitle("coloured by age and day")

################################################################################
#Coloured by individual
ggplot(data = coord_PCA_tb, aes(x=pc1, y=pc2)) +
  geom_point(aes(color = indiv), size=1) +
  xlab(coord_PCA_calc[["labs1"]]) +
  ylab(coord_PCA_calc[["labs2"]]) +
  theme_bw() #+
  #ggtitle("coloured by individual")

#coord_cartesian(xlim = c(,), ylim = c(,))
#