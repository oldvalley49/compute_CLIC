library(Signac)
library(Seurat)
library(GenomicRanges)
library(glue)
library(EnsDb.Hsapiens.v86)
library(rtracklayer)

### HUMAN

# we generate seurat objects for each multiomics sample, filtering cells according to their qc metric

# argument after script is multiome sample accession
args = commandArgs(trailingOnly = TRUE)
sample_multiome_acc = args[1]

# annotations UCSC style
annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

# create Seurat Object and filter based on qc metrics

inputdata.10x = Read10X_h5(glue('data/human/{sample_multiome_acc}/filtered_feature_bc_matrix.h5'))
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# RNA assay

sample_obj = CreateSeuratObject(counts = rna_counts, assay = "RNA", project = sample_multiome_acc)
sample_obj[["percent.mt"]] <- PercentageFeatureSet(sample_obj, pattern = "^MT-")

# ATAC assay

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

sample_obj[["ATAC"]] = CreateChromatinAssay(counts = atac_counts,
                                            fragments = glue('data/human/{sample_multiome_acc}/atac_fragments.tsv.gz'),
                                            sep = c(":", "-"),
                                            annotation = annotations)
DefaultAssay(sample_obj) = "ATAC"

sample_obj = NucleosomeSignal(object = sample_obj, verbose = F)
sample_obj = TSSEnrichment(object = sample_obj, verbose = T, fast = T)  

# compute gene activity using Signac
gene.activities <- GeneActivity(sample_obj, features = rownames(rna_counts))
sample_obj[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# filter by QC metrics

DefaultAssay(sample_obj) <- "RNA"
sample_obj = subset(
	x = sample_obj, 
	subset = nCount_ATAC < 60000 &
    nCount_ATAC > 500 &
	nCount_RNA < 40000 &
	nCount_RNA > 500 &
	nucleosome_signal < 2 &
	TSS.enrichment > 1 &
	percent.mt < 10
)

# save remapped Seurat Object
saveRDS(sample_obj, glue("output/init_obj/human/{sample_multiome_acc}.rds"))


