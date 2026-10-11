#!/bin/bash
#SBATCH --job-name=hfd_cassette_B
#SBATCH --partition=general
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/cassetteB_%A_%a.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/cassetteB_%A_%a.err

# Step 2.4 stage B: place reads on the reconstructed cassettes.
#
# Two streams per sample, kept separate:
#   aligned  - reads currently assigned to the four v1 transgene contigs
#   unmapped - every unplaced-unmapped read (region '*', so the 29 GB BAM is
#              never decompressed end to end)
#
# The unmapped stream is the one that matters. 2_2 already scanned it, but only
# for k-mers of the four v1 contigs, so any read from the promoter, the intron,
# the ORF-to-WPRE linker, the ITRs, or the parts of WPRE/SV40pA that differ
# between constructs was invisible to it. Here the whole ITR-to-ITR cassette is
# the search space.

set -euo pipefail

module load samtools/1.18-gcc-12.3.0

PROJ=/scratch/user/uqschauq/Insular
MANIFEST=$PROJ/Data/sample_metadata.csv
COUNTS=$PROJ/Data/Counts
BLOCKS=$PROJ/Results/2_4_cassette_blocks
OUTDIR=$PROJ/Results/2_4_cassette_coverage
SCANNER=$PROJ/Code/2_4_cassette_coverage.py

mkdir -p "$OUTDIR"

for f in cassettes.fa cassette_basemask.tsv; do
    [[ -s $BLOCKS/$f ]] || { echo "FATAL: stage A output $BLOCKS/$f missing"; exit 1; }
done

# fastq_prefix of rows with fastq_ok==True; the manifest has CRLF line endings
SAMPLE=$(tr -d '\r' < "$MANIFEST" | awk -F, -v n="$SLURM_ARRAY_TASK_ID" '
  NR==1 { for (i=1; i<=NF; i++) { if ($i=="fastq_prefix") p=i; if ($i=="fastq_ok") k=i } ; next }
  $k=="True" { c++; if (c==n) { print $p; exit } }')

if [[ -z "$SAMPLE" ]]; then
  echo "No sample for array index $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

BAM=$COUNTS/$SAMPLE/outs/possorted_genome_bam.bam
[[ -s $BAM ]] || { echo "FATAL: missing $BAM"; exit 1; }

echo "[$(date)] $SAMPLE: stage B on $BAM"

echo "[$(date)] $SAMPLE: aligned-to-transgene stream"
samtools view -@ 4 "$BAM" EGFP tdTomato EBFP2 AAV_WPRE_SV40pA_shared \
  | python3 "$SCANNER" "$SAMPLE" aligned "$OUTDIR"

echo "[$(date)] $SAMPLE: unmapped stream"
samtools view -@ 4 "$BAM" '*' \
  | python3 "$SCANNER" "$SAMPLE" unmapped "$OUTDIR"

echo "[$(date)] $SAMPLE done"
cat "$OUTDIR/${SAMPLE}_aligned_summary.tsv" "$OUTDIR/${SAMPLE}_unmapped_summary.tsv"
