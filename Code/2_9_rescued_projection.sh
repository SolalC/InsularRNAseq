#!/bin/bash
#SBATCH --job-name=hfd_proj29
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/proj29_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/proj29_%j.err

# Step 2.9: fold the 2.8 rescued molecules into the per-cell transgene counts and
# re-run the projection classification, before vs after, per animal.

set -euo pipefail
python3 /scratch/user/uqschauq/Insular/Code/2_9_rescued_projection.py
