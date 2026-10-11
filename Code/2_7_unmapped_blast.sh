#!/bin/bash
#SBATCH --job-name=hfd_ublast
#SBATCH --partition=general
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/ublast_%A_%a.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/ublast_%A_%a.err

# Step 2.7: BLAST every unmapped read against the three full transgene cassettes.
#
# Why this is not a repeat of 2_2 / 2_4. Both of those scanned for EXACT 31-mers:
# a read only counted if it carried a perfect 31 bp window. That misses
#   - reads with scattered sequencing errors (a mismatch every <31 bp),
#   - reads that only partly overlap the cassette (integration junctions,
#     cassette-to-genome chimeras, adapter/polyA runoff),
#   - construct sequence that diverges from the plasmid map (recombination,
#     ITR rearrangement, a variant WPRE),
#   - anything under 31 clean bases.
# blastn with word_size 11 finds all of those, so this is a genuine sensitivity
# floor for "is there transgene sequence among the reads that mapped nowhere".
#
# Scope caveat: this pool is reads that failed to map to mouse AND to the v2
# transgene contigs. Transgene reads that mapped correctly are not in it, so
# this bounds what was MISSED, it does not recount what was found.
#
# Database is the full ITR-to-ITR cassettes from 2_4 (unmasked, 10,079 bp over
# three sequences), not the v2 reference contigs - the point is to search the
# complete construct, including sequence that is masked to N in the reference.

set -euo pipefail

module purge
module load samtools/1.18-gcc-12.3.0
module load blast+/2.14.1-gompi-2023a

PROJ=/scratch/user/uqschauq/Insular
COUNTS=$PROJ/Data/Counts_v2
MANIFEST=$PROJ/Data/sample_metadata.csv
CASSETTES=$PROJ/Results/2_4_cassette_blocks/cassettes.fa
OUT=$PROJ/Results/2_7_unmapped_blast
DBDIR=$OUT/db
CORES=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "$OUT" "$DBDIR"

# sample_metadata.csv is CRLF; strip \r and look the column up by header name.
SAMPLE=$(awk -F',' -v n="$SLURM_ARRAY_TASK_ID" '
  NR==1 { for (i=1;i<=NF;i++) { h=$i; gsub(/\r/,"",h); if (h=="fastq_prefix") c=i }; next }
  NR==n+1 { v=$c; gsub(/\r/,"",v); print v }' "$MANIFEST")

[[ -n $SAMPLE ]] || { echo "FATAL: no sample for array index $SLURM_ARRAY_TASK_ID"; exit 1; }

BAM=$COUNTS/$SAMPLE/outs/possorted_genome_bam.bam
[[ -s $BAM     ]] || { echo "FATAL: no BAM at $BAM"; exit 1; }
[[ -s $BAM.bai ]] || { echo "FATAL: no BAM index at $BAM.bai"; exit 1; }

# One shared BLAST database. Task 1 builds it; the others wait for it rather
# than racing to write the same files.
DB=$DBDIR/cassettes
if [[ $SLURM_ARRAY_TASK_ID -eq 1 ]]; then
  makeblastdb -in "$CASSETTES" -dbtype nucl -out "$DB" -title cassettes >"$DBDIR/makeblastdb.log"
  touch "$DBDIR/.ready"
else
  for _ in $(seq 1 120); do [[ -e $DBDIR/.ready ]] && break; sleep 5; done
  [[ -e $DBDIR/.ready ]] || { echo "FATAL: BLAST db not ready after 10 min"; exit 1; }
fi

echo "[$(date)] $SAMPLE: blasting unmapped reads (${CORES} threads)"

# Unplaced-unmapped reads are fetched through the BAM index (region '*'), so the
# ~29 GB BAM is never decompressed end to end. NR makes the query id unique even
# if a read name recurs. dust is left on: low-complexity polyA/polyG runoff would
# otherwise carry the hit table.
samtools view -@ 4 "$BAM" '*' \
  | awk '{ print ">r" NR "_" $1; print $10 }' \
  | blastn \
      -db "$DB" \
      -task blastn \
      -word_size 11 \
      -evalue 1e-5 \
      -max_target_seqs 5 \
      -num_threads "$CORES" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  > "$OUT/${SAMPLE}_hits.tsv"

echo "[$(date)] $SAMPLE: $(wc -l < "$OUT/${SAMPLE}_hits.tsv") HSPs"
echo "[$(date)] done"
