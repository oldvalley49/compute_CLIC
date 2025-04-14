library(Seurat)
library(Signac)
library(glue)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(stringr)
options(future.globals.maxSize = 1e9)
library(tidyverse)
library(stringr)
library(harmony)
# we integrate samples originating from the same tissue using harmony


# obtain samples and corresponding tissues
all_samples_tissue_df = read.csv("metadata/ENCODE_metadata_mouse.csv")

pdf("plots/mouse_harmonyIntegrated.pdf", width = 14, height = 6)

# get list of tissues
tissue_list = unique(all_samples_tissue_df$tissue_source)
cluster_num <- list()
# iterate through tissues
for (tissue in tissue_list) {

	# obtain samples pertaining to same tissue
	samples_per_tissue_list = all_samples_tissue_df$multiomics_series[all_samples_tissue_df$tissue_source == tissue]
	
	# obtain samples
	samples_per_tissue_filepaths = lapply(samples_per_tissue_list, function(x) glue("output/init_obj/mouse/{x}.rds"))

	tissue_multiomics_obj_list = lapply(samples_per_tissue_filepaths, readRDS)

	tissue_merged_obj = merge(x = tissue_multiomics_obj_list[[1]], y = tissue_multiomics_obj_list[-1], add.cell.ids = samples_per_tissue_list)

	# normalization
	tissue_merged_obj <- NormalizeData(object = tissue_merged_obj)
	tissue_merged_obj <- FindVariableFeatures(object = tissue_merged_obj)
	tissue_merged_obj <- ScaleData(object = tissue_merged_obj)
	tissue_merged_obj <- RunPCA(tissue_merged_obj)

	#clustering on unintegrated samples
	# tissue_merged_obj <- FindNeighbors(tissue_merged_obj, reduction = "pca", dims = 1:30)
	# tissue_merged_obj <- FindClusters(tissue_merged_obj, resolution = 0.5, cluster.name = glue("Unintegrated_clusters"))

	# unintegrated umap
	tissue_merged_obj <- RunUMAP(tissue_merged_obj, reduction = "pca", dims = 1:30, reduction.name = 'umap.unintegrated')

	# harmony integration
	tissue_merged_obj <- RunHarmony(tissue_merged_obj, group.by.vars='orig.ident')
	
	# clustering on harmony integrated samples
	tissue_merged_obj <- FindNeighbors(tissue_merged_obj, reduction = "harmony", dims = 1:30)
	tissue_merged_obj <- FindClusters(tissue_merged_obj, resolution = 0.5, cluster.name = "Harmony_clusters")

	# save info about number of cells in each cluster
	clusters_list <- tissue_merged_obj[["Harmony_clusters"]]
	clusters_counts <- table(clusters_list)
	cluster_num[[tissue]] <- clusters_counts

	# harmony integrated umap
	tissue_merged_obj <- RunUMAP(tissue_merged_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
	
	# plot original unintegrated umap by sample
	
	p1 = DimPlot(tissue_merged_obj,
             	reduction = "umap.unintegrated",
             	group.by = c("orig.ident"),
             	combine = TRUE, label.size = 2
	) + ggtitle(glue("Unintegrated {tissue} samples"))

	# plot harmony integrated umap by sample
	p2 <- DimPlot(tissue_merged_obj,
                reduction = "umap.harmony",
                group.by = c("orig.ident"),
                combine = TRUE, label.size = 2
	) + ggtitle(glue("Harmony Integrated {tissue} samples"))

	p3 = p1 + p2
	
	print(p3)

	# save merged integrated tissue sample
	tissue_underscore = str_replace_all(tissue, " ", "_")
	saveRDS(tissue_merged_obj,glue("output/tissue_harmonized/mouse/{tissue_underscore}.rds"))
}	

saveRDS(cluster_num, glue("output/tissue_harmonized/cluster_info_mouse.rds"))

dev.off()


