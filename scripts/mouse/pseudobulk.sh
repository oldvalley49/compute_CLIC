#!/bin/bash

#SBATCH --job-name=mouse_pseudobulk
#SBATCH --time=3:00:00
#SBATCH --mem=150G

module load R

Rscript --no-save --no-restore src/mouse/pseudobulk.R














