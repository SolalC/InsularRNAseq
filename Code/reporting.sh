#!/bin/bash
#SBATCH --job-name=hfd_reporting
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=300G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/reporting_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/reporting_%j.err

# Figures for the write-up. Housekeeping, not a pipeline step -- see the header
# of reporting.R. Nothing downstream depends on it and it writes no object.
#
#   sbatch reporting.sh                          # both sections, default genes
#   sbatch reporting.sh doublet                  # violins only (~120G is enough)
#   sbatch reporting.sh decon Plp1,Mbp,Snap25    # pre/post feature plots only
#   sbatch reporting.sh decon Plp1,Mbp /scratch/user/uqschauq/Insular/Data/Seurat/HFDseq_v2_umap.rds
#
# The memory here is sized for the "decon" section, which holds the 17 GB UMAP
# object and re-reads a 1.8M-barcode raw matrix at the same time.

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none} args=${*:-none}"
Rscript --vanilla "$PROJ/Code/reporting.R" "${1:-all}" "${2:-}" "${3:-}"
echo "[$(date)] done"
