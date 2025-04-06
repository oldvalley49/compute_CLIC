#!/bin/bash

#SBATCH --job-name=init_seurat_human
#SBATCH --output=logs/init_seurat_human/r_method_%A_%a.out
#SBATCH --error=logs/init_seurat_human/r_method_%A_%a.err
#SBATCH --array=11
#SBATCH --time=5:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1


# NOTE: jobs 8--15 needs  time (tried with 3 hrs)

module load R

# read parameters from human ENCODE metadata
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/ENCODE_human.txt)

# Parse parameters
read ENCODE_ACCESSION <<< "$PARAMS"

echo "ENCODE Accession: $ENCODE_ACCESSION"

Rscript --no-save --no-restore src/human/init_seurat.R $ENCODE_ACCESSION

