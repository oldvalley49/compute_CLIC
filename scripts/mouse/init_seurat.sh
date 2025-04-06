#!/bin/bash

#SBATCH --job-name=init_seurat_mouse
#SBATCH --output=logs/init_seurat_mouse/r_method_%A_%a.out
#SBATCH --error=logs/init_seurat_mouse/r_method_%A_%a.err
#SBATCH --array=1-40
#SBATCH --time=1:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1

module load R

# read parameters from mouse ENCODE metadata
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/ENCODE_mouse.txt)

# Parse parameters
read ENCODE_ACCESSION <<< "$PARAMS"

echo "ENCODE Accession: $ENCODE_ACCESSION"

Rscript --no-save --no-restore src/mouse/init_seurat.R $ENCODE_ACCESSION

