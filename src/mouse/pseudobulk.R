library(Seurat)
library(Signac)
library(glue)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(stringr)
options(future.globals.maxSize = 1e9)

# obtain samples and corresponding tissues
all_samples_tissue_df = read.csv("metadata/ENCODE_metadata_mouse.csv")

# get list of tissues
tissue_list = unique(all_samples_tissue_df$tissue_source)

for (tissue in tissue_list) {

	tissue_underscore = str_replace_all(tissue, " ", "_")

	# load harmony integrated tissue sample with cluster info based on harmony embedding
	tissue_samples_clustered = readRDS(glue("output/tissue_harmonized/mouse/{tissue_underscore}.rds"))
	
	# pseudobulk counts by clusters
	tissue_samples_pseudobulked = AggregateExpression(tissue_samples_clustered, group.by=c("Harmony_clusters"))
	
	# save pseudobulked tissue samples in new file
	saveRDS(tissue_samples_pseudobulked, file = glue("output/harmony_pseudobulk/mouse/{tissue_underscore}.rds"))
}

