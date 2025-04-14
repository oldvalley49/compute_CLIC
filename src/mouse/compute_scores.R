library(dplyr)
library(preprocessCore)
library(Matrix)
library(ggplot2)
library(glue)
### read all the pseudobulked data

# obtain samples and corresponding tissues
all_samples_tissue_df = read.csv("metadata/ENCODE_metadata_mouse.csv")

# get list of tissues
tissue_list = unique(all_samples_tissue_df$tissue_source)

data_list <- c()

for (tissue in tissues_list) {
    data_list <- c(data_list, readRDS(glue("output/harmony_pseudobulk/mouse/{tissue}.rds")))
}

# combined data


