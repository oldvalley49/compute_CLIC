#!/bin/bash

#SBATCH --job-name=human_harmonize
#SBATCH --time=1-00:00:00
#SBATCH --mem=200G

module load R

Rscript --no-save --no-restore src/human/integrate_harmony.R














