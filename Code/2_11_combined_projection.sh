#!/bin/bash
#SBATCH --job-name=hfd_proj211
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/proj211_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/proj211_%j.err

# Step 2.11: projection calls from the v2 matrix plus both rescue tracks,
# reported as base -> +BLAST rescue -> +antisense so either can be dropped.

set -euo pipefail
python3 /scratch/user/uqschauq/Insular/Code/2_11_combined_projection.py
