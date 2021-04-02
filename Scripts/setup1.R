#Clear
#...variables
rm(list=ls())
#...console
cat("\014\n")
#...graphs
tryCatch(dev.off(), error = function(e) {NULL})
dev.new()

library(tidyverse)
library(MAST)

################################################################################

#Data sub-folder inside the project
raw <- readr::read_csv("./Data/GEO_supporting_processed_data_file_raw_count_matrix.csv")