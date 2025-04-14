library(rtracklayer)
library(GenomicRanges)
library(Matrix)
library(ggplot2)
library(dplyr)
library(preprocessCore)
library(glue)
library(stringr)
library(Seurat)
all_samples_tissue_df = read.csv("metadata/ENCODE_metadata_human.csv")

# get list of tissues
tissue_list = unique(all_samples_tissue_df$tissue_source)
tissue_list <- lapply(tissue_list, function(x) str_replace_all(x, " ", "_"))

# load pseudobulked data

tissue_pseudobulked_list <- list()
for (tissue in tissue_list) {
    data <- readRDS(glue("output/harmony_pseudobulk/human/{tissue}.rds"))
    data <- lapply(data, function(x) {
        colnames(x) <- paste0(tissue, "_", colnames(x))
        return(x)
        })
    tissue_pseudobulked_list[[tissue]] <- data
}

rna_pseudobulk_list <- lapply(tissue_pseudobulked_list, function(x) x[["RNA"]])
activity_pseudobulk_list <- lapply(tissue_pseudobulked_list, function(x) x[["ACTIVITY"]])

# combine the data into one matrix
rna_counts <- do.call(cbind, rna_pseudobulk_list)
activity_counts <- do.call(cbind, activity_pseudobulk_list)

### Normalize Count Matrices

# subset by genes that have both expression and gene activity
rna_counts <- rna_counts[rownames(activity_counts), ]

rna_counts <- NormalizeData(rna_counts)
rna_data <- as.data.frame(normalize.quantiles(as.matrix(rna_counts), keep.names = TRUE))

activity_counts <- NormalizeData(activity_counts)
activity_data <- as.data.frame(normalize.quantiles(as.matrix(activity_counts), keep.names = TRUE))

saveRDS(rna_data, "output/normalized_data/human_rna.rds")
saveRDS(activity_data, "output/normalized_data/human_activity.rds")

rna_data <- as.matrix(rna_data)
activity_data <- as.matrix(activity_data)
genes <- rownames(rna_data)
correlation <- data.frame(matrix(NA, nrow = length(genes), ncol = 1 ))
colnames(correlation) <- c("pearson_correlation")
rownames(correlation) <- genes
for (gene in genes){
    correlation[gene, "pearson_correlation"] <- cor(rna_data[gene, ], activity_data[gene, ], method = "pearson")
}

write.csv(correlation, "output/scores/human_scores.csv")




