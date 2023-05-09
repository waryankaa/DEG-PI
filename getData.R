#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(tidyverse)

rm(list = ls())

name <- read.delim("/Users/abbywaryanka/seniorComp/Streamlit/keyword.txt", header = FALSE)
# load series and platform data from GEO
gset <- getGEO(name, GSEMatrix =TRUE, AnnotGPL=TRUE)
ending <- "_series_matrix.txt.gz"
type <- paste(name,ending)
type <- str_replace_all(string=type, pattern=" ", repl="")
platform <- gset[[type]]@annotation

if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# categorize samples
sample_type_table <- select(gset@phenoData@data,title,geo_accession)
write.csv(sample_type_table, "/Users/abbywaryanka/seniorComp/Streamlit/sample_type.csv", row.names=FALSE)







