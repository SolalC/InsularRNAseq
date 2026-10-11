#!/bin/bash
#SBATCH --job-name=hfd_multiorf
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/multiorf_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/multiorf_%j.err

# Step 2.12: how many nuclei carry more than one reporter, which region pairs,
# and is that above the coincidence expectation for multiplets/ambient?

set -euo pipefail
python3 /scratch/user/uqschauq/Insular/Code/2_12_multiorf_analysis.py
