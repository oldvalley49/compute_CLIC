#!/bin/bash

#SBATCH --job-name=mouse_harmonize
#SBATCH --time=10:00:00
#SBATCH --mem=150G

module load R

Rscript --no-save --no-restore src/mouse/integrate_harmony.R














