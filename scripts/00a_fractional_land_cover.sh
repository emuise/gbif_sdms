#!/bin/bash
#SBATCH --account=def-reyna
#SBATCH --array=1-19
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03S:00:00
#SBATCH --job-name=lc_frac
#SBATCH --output=%x-%A_%a.out

cd "$SLURM_SUBMIT_DIR"
module load StdEnv gcc r gdal geos proj
Rscript scripts/00a_fractional_land_cover.R