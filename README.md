# CLIC_score
Computing correlation score for scRNA-seq scATAC-seq integration using ENCODE data


TODO: 
1. test init, qc, harmonizing, and pseudobulking -> generate data, check plots for integration
2. normalize pseudobulk data (must be normalized by gene width -> look at GeneActivity() code in Signac) -> quantile normalize
3. generate plots for integration score to use in manuscript

# Reproducing

1. Download the ENCODE raw sequencing data
2. Run through Cell Ranger and output in `data/{species}/` for human and mouse, respectively
3. For each species, run 
