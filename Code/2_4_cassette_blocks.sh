#!/bin/bash
#SBATCH --job-name=hfd_cassette_A
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/cassetteA_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/cassetteA_%j.err

# Step 2.4 stage A: reconstruct the packaged AAV cassette (ITR to ITR) for each
# of the three constructs from p24/p131/p140, partition the three cassettes by
# sequence sharing, and measure how much of each cassette the v1 reference left
# unmodelled. Produces the k-mer index and per-base mask used by stage B.
#
# Uses the real plasmid maps, NOT Data/Plasmid/plasmidConstruct.FA - the
# KJ081243.1 record in that file is only 25.8% nucleotide-identical to p131's
# codon-diversified tdTomato.

set -euo pipefail

PROJ=/scratch/user/uqschauq/Insular
OUTDIR=$PROJ/Results/2_4_cassette_blocks

mkdir -p "$OUTDIR" "$PROJ/Code/logs"

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"

for f in p24.gb p131.gb p140.gb; do
    [[ -s $PROJ/Data/Plasmid/$f ]] || { echo "FATAL: missing Data/Plasmid/$f"; exit 1; }
done
[[ -s $PROJ/Data/Reference/transgenes/transgenes.fa ]] \
    || { echo "FATAL: missing the v1 transgene contigs"; exit 1; }

python3 "$PROJ/Code/2_4_cassette_blocks.py"

echo "[$(date)] stage A done"
ls -la "$OUTDIR"
